## Python Code for control of plasmonic dPCR device V1.0
## Christian D. Ahrberg
## Sogang University - November 2019

#--------------------------------------------------------------------------

# Importing required packages
import time
import os
import RPi.GPIO as GPIO
import tkinter as tk
import sqlite3
import datetime
from scipy.signal import butter, lfilter
import numpy as np
from numpy import convolve

# -------------------------------------------------------------------------

# Function for reading value from ADC

# Create function for measurment from ADC (8 channels 0 to 7)
def readadc(adcnum, clockpin, mosipin, misopin, cspin):
    # Checking for correct channel
    if((adcnum > 7) or (adcnum < 0)):
        return -1
    GPIO.output(cspin,True)

    GPIO.output(clockpin, False)    # Start clock low
    GPIO.output(cspin,False)        # bring CS low

    commandout = adcnum
    commandout |= 0x11              # Check, startbit + single end bit
    commandout <<= 3                # Only sending 5 bits
    for i in range(5):
        if(commandout & 0x80):
            GPIO.output(mosipin,True)
        else:
            GPIO.output(mosipin,False)
        commandout <<= 1
        GPIO.output(clockpin, True)
        GPIO.output(clockpin, False)

    adcout = 0
    # read one empty bit, one null bit and 10 ADC bits
    for i in range(15):
        GPIO.output(clockpin, True)
        GPIO.output(clockpin, False)
        adcout <<= 1
        if (GPIO.input(misopin)):
            adcout |= 0x1

    GPIO.output(cspin, True)

    adcout >>= 1            # First bit is null, so drop it
    adcout = adcout - 2*4095
    return adcout
# --------------------------------------------------------------------------

# Function for converting ADC to voltage

def ADCtoV(ADC, Vref):
    
    # Function converting the ADC reading of MCP3008 to a voltage

    # Definition of inputs:
    # ADC   - reading form MCP3008
    # Vref  - Reference voltage to MCP3008 [V]
    
    ADC = float(ADC)    # making sure ADC reading is a float\
    Voltage = ADC/4095 * Vref   # conversion to voltage
    return Voltage
# --------------------------------------------------------------------------

# PID controller function

def PIDcont(Setpoint, MeasPoint,Kp,Ki,Kd,t,i0,e0):

    # PID controller function returns output value of PID controller (uPID), error of current itteration (e)
    # , and the current value of the integral for the Integral controller (I)

    # make sure function is called every t seconds

    # Definittion of inputs:
    # Setpoint = Desired value, here temperature
    # MeasPoint = Measured value, here temperature
    # Kp = Value for proptoional controller
    # Ki = Value for integral controller
    # Kd = Value for differential controller
    # t = Sampling time [s]
    # i0 = initial state of integrator (default =0)
    # e0 = Initial error (default = 0)

    # Initialisation
    Ii_prev = float(i0)     # Previous integration
    e_prev = float(e0)      # Previous error

    Kp = float(Kp)
    Ki = float(Ki)
    Kd = float(Kd)
    t = float(t)

    # Actual contoller

    
    e = Setpoint - MeasPoint    # Callculating current error

    # Proportional:
    uP = Kp * e                 # Proportional part

    # Integral:
    I = Ii_prev + t * e         # Integral
    uI = Ki * I                 # Integral part
    Ii_prev = I                 # Updating integral for next interation

    # Differential
    dedt = (e - e_prev) / t     # Differential of error
    uD = Kd * dedt              # Differential part
    e_prev = e                  # Updating error

    # Overall controller output

    uPID = uP + uI + uD         # Overall PID output

    return (uPID,e,I)

# --------------------------------------------------------------------------

# Function for converting voltage to temperature

def VtoT(Vin, Vsupp, R1, R2, R3, a, b):

    # Function converting voltage from a Wheatstone bridge to temperature measurement

    # Definition of inputs:
    # Vin   - Measured voltage of Wheatstone bridge [V]
    # Vsupp - Voltage supplied to Wheatstone bridge [V]
    # R1    - First resistance [Ohm]
    # R2    - Adjustable resistanc in parralel with measured resistance [Ohm]
    # R3    - Resistance in parallel with R1 [V]
    # a     - Factor for converting Resistance to T (T = (R-a)/b [Ohm]
    # b     - Factor for converting Resistance to T (T = (R-a)/b [Ohm/C]

    alpha = R2 / (R1 + R2) - Vin/Vsupp  # Dividing callculation into two blocks
    Rmes = alpha * R1/(1 - alpha)       # Callculating the resistance of temperature sensor
    T = (Rmes - a) / b                  # Converting resistence to temperatrue according to callibration
        
    return T

#---------------------------------------------------------------------------

# main function with all the code for cycling

def writeTprof(Ncyc,tHS,Ths,TDE,Tde,tAN,Tan,tEX,Tex):

    # Initialsing, setting up GPIO labelling
    GPIO.setmode(GPIO.BCM)
    DEBUG = 1
    
    # SQL operations with database
    conn = sqlite3.connect('CyclingData.db')    # Connecting to database
    c = conn.cursor()                           # Creating cursor for access
    # Creating table with experiments index (if it does not exist)
    c.execute('CREATE TABLE IF NOT EXISTS experiments (Date varchar(255), Dataname varchar(255))')
    # Creating table for Temperature profile
    c.execute('CREATE TABLE IF NOT EXISTS tempprofile (Timename varchar(255), Time varchar(255), Tempname varchar(255),Temp varchar(255))')
    conn.commit()

    # Creating table in SQL database to save data to
    date = datetime.date.today()                        # Determining current data
    rows = c.execute('''select * from experiments''')   # Finding previous number of experiments
    rows = rows.fetchall()
    rows = len(rows)+1                                  # Incrementing by one to create new unique index
    dataname = 'data_'+str(rows)                        # Unique name for table containging cycling data
    c.execute('Insert INTO experiments (Date, Dataname) VALUES (?,?)',(date,dataname))  # Adding thermal data table to db
    conn.commit()                                       # Commiting change to db

    # Creating table for Experimental data
    c.execute('CREATE TABLE IF NOT EXISTS ' + dataname + ' (Time, ADC, Volt, MeasT, SetT, PID, PWM, Temp)')

    # Importing experimental constants from the database
    variables = c.execute('SELECT Value FROM Variables')
    variables = [float(item[0]) for item in variables.fetchall()]
    # Constants as imported from database
    Freq = variables[0]     # Number of measurements per second [Hz]
    Vref = variables[1]     # Reference voltage from ADC [V]
    Vsupp = variables[2]    # Voltage supplied to wheatsonte bridge [V]
    R1 = variables[3]       # Resistance for wheatstone bridge [Ohm]
    R2 = variables[4]       # Adjustable resistance for Wheatstone bridge [Ohm]
    R3 = variables[5]       # Resistance for wheatstone bridge [Ohm]
    a = variables[6]        # Factor for resistance to temperature conversion [Ohm]
    b = variables[7]        # Factor for resistance to temperature conversion [Ohm/C]
    Kp = variables[8]       # Value for K part
    Ki = variables[9]       # Value for I part
    Kd = variables[10]      # Value for d part
    fc = variables[11]      # Cut of frequency in Hz for low pass filter
    freq_PWM = variables[12]# Frequency for PWM in Hz
    maxPID = variables[13]  # Value for conversion of PID to PWM
    fan_threshold = variables[14] # threshold for PID value, below which fan will turn on

    # Initialising parameters for cycling
    e = 0           # Inital error for PID
    I = 0           # Inital integral for PID
    j = 0           # Counter for cycles
    i = 0           # Counter for steps
    intTime = 0     # Internal timer for step times
    setTime = 0     # Target time for step
    Time1 = 0       # Overall elappsed time
    fanon = 0       # is the cooling fan on? 0 for no, 1 for yes
    dcycle_PWM = 50 # Default duty cycle for PWM in %
    fs = float(Freq)    # Sampling frequency for low pass filter
    w = fc / (fs / 2)   # Normalize the frequency for low pass filter

    # Making sure everything is a float
    Vref  = float(Vref)
    Vsupp = float(Vsupp)
    R1    = float(R1)
    R2    = float(R2)
    R3    = float(R3)
    a     = float(a)
    b     = float(b)
    Freq  = float(Freq)

    # Defingin pin connections
    SPICLK  = 18     # For ADC
    SPIMISO = 23     # For ADC
    SPIMOSI = 24     # For ADC
    SPICS   = 25     # For ADC
    PINPWM  = 27     # Pin for PWM
    PINFAN  = 17     # Pin for active cooling with fan

    # Define Pin connection to ADC 3008
    potentiometer_adc = 1 #Ch2

    # Writing values of used temperature profile to Tprof in database
    c.execute('UPDATE tempprofile SET Time = ? WHERE Timename = ?',(Ncyc.get(),'Cycle'))
    c.execute('UPDATE tempprofile SET Time = ? WHERE Timename = ?',(tHS.get(),'Hot Start'))
    c.execute('UPDATE tempprofile SET Temp = ? WHERE Timename = ?',(Ths.get(),'Hot Start'))
    c.execute('UPDATE tempprofile SET Time = ? WHERE Timename = ?',(TDE.get(),'Denaturation'))
    c.execute('UPDATE tempprofile SET Temp = ? WHERE Timename = ?',(Tde.get(),'Denaturation'))
    c.execute('UPDATE tempprofile SET Time = ? WHERE Timename = ?',(tAN.get(),'Annealing'))
    c.execute('UPDATE tempprofile SET Temp = ? WHERE Timename = ?',(Tan.get(),'Annealing'))
    c.execute('UPDATE tempprofile SET Time = ? WHERE Timename = ?',(tEX.get(),'Extension'))
    c.execute('UPDATE tempprofile SET Temp = ? WHERE Timename = ?',(Tex.get(),'Extension'))
    conn.commit()
    
    # Setup SPI interface pins
    GPIO.setup(SPIMOSI, GPIO.OUT)
    GPIO.setup(SPIMISO, GPIO.IN)
    GPIO.setup(SPICLK, GPIO.OUT)
    GPIO.setup(SPICS, GPIO.OUT)
    GPIO.setup(PINPWM, GPIO.OUT)
    GPIO.setup(PINFAN, GPIO.OUT)

    
    # Extracting values from imported text file
    profile_time = c.execute('SELECT Time FROM tempprofile')
    profile_time = [int(item[0]) for item in profile_time.fetchall()]
    profile_temp = c.execute('SELECT Temp FROM tempprofile')
    profile_temp = [int(item[0]) for item in profile_temp.fetchall()]

    tDE = profile_time[0] # Time for denaturation
    Tde = profile_temp[0] # Temperature for denaturation step       
    tAN = profile_time[1] # Time for annealing step
    Tan = profile_temp[1] # Temperature for annealing step
    tEX = profile_time[2] # Time for extension step
    Tex = profile_temp[2] # Temperature for extesnion step
    Nc  = profile_time[3] # Number of cycles
    tHS = profile_time[4] # Time for hot start
    Ths = profile_temp[4] # Temperatrue for hot start

    # Calculating number of measurments and waitung time
    Waitt = 1/Freq # Waiting time between two measurments [s]

    # Writing all used varibales to database for later reference
    c.execute('''CREATE TABLE IF NOT EXISTS ExperimentsVariables
              (Date, Experimentname, Freq, Vref, Vsupp, R1, R2, R3, a, b, Kp, Ki, Kd, fc, freq_PWM, maxPID, fan_threshold, 
              tDE, Tde1,tAN, Tan1, tEX, Tex1, Nc, tHS, Ths1)''')
    c.execute('''INSERT INTO  ExperimentsVariables VALUES
             (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,?, ?, ?, ?, ?, ?, ?, ?, ?
             )''',(str(date), str(dataname), Freq, Vref, Vsupp, R1, R2, R3, a, b, Kp, Ki, Kd, fc, freq_PWM, maxPID, fan_threshold, tDE, Tde, tAN, Tan, tEX, Tex, Nc, tHS, Ths))
    conn.commit()

    # Setting up PWM
    func_PWM = GPIO.PWM(PINPWM, freq_PWM)

    # Measurment loop

    func_PWM.start(0) # Starting PWM

    while j <= Nc:

        if (j == 0 & i != 66.6): # Hot start, random number so true for all i, if j = 0
            setT = Ths
            setTime = tHS
            i = 2 # So once timer is up this gets updated to cycling
        elif (j >= 0 and i == 0):     # Cycling: Denaturation
            setT = Tde
            setTime = tDE
        elif (j >= 0 and i == 1):     # Cycling: Annealing
            setT = Tan
            setTime = tAN
        elif (j >= 0 and i == 2):     # Cycling: Extension
            setT = Tex
            setTime = tEX
        
        # Reading ADC value
        adcout = readadc(potentiometer_adc, SPICLK, SPIMOSI, SPIMISO, SPICS)
        # Converting to voltage
        voltage = ADCtoV(adcout,Vref)
        # Converting to temperature
        Temperature = VtoT(voltage, Vsupp, R1, R2, R3, a, b)
        # Determining time
        Time1 = Time1 + Waitt
        # Using PID
        (uPID,e,I) = PIDcont(setT, Temperature,Kp,Ki,Kd,Waitt,e,I)

        # Determining new Duty cycle to achieve heating / cooling
        dcycle_PWM = uPID / maxPID * 100 # Calculating new dutycycle according to PID
        # Making sure new Duty cycle is maximum 100% and minimum 0
        if dcycle_PWM > 100:
            dcycle_PWM = 100
        elif dcycle_PWM < 0:
            dcycle_PWM = 0

        # Updating Duty Cycle
        func_PWM.ChangeDutyCycle(dcycle_PWM)
        
        # Turn cooling fan on if PID << 0
        # Turn off if PID > 0
        if (uPID <= fan_threshold) and (fanon == 0): # turn on fan
            fanon = 1
            GPIO.output(PINFAN,True)
        elif (uPID <= fan_threshold) and (fanon == 1): # turn on fan
            pass
        else:  # turn off fan
            fanon = 0
            GPIO.output(PINFAN,False) 
           
        # Writing data to database
        c.execute('INSERT INTO ' +dataname + ' VALUES ('+ str(Time1) +',' + str(adcout) +', ' + str(voltage) + ', ' + str(Temperature) + ', ' + str(setT) + ', ' + str(uPID) + ', ' + str(dcycle_PWM) + ', 0)')
        
        time.sleep(Waitt)

        # Updating timers and moving to next step / cycle if required
        intTime = intTime + Waitt

        if intTime >= setTime:  # Time for step ellapsed
            i = i+1             # Moving on to next step
            intTime = 0         # Reseting time
            e = 0               # Reseting error for PID
            I = 0               # Reseting integral for PID
            
            if i >= 3:          # Moving to next cycle
                i = 0           # Reseting step counter
                j = j+1         # Moving to next cycle
                e = 0           # Reseting error for PID
                I = 0           # Reseting integral for PID
                
    # Clearing GPIOs and stopping LEDs after all cycles
    func_PWM.stop() # Stoping PWM
    GPIO.cleanup()

    # Filtering temperature data to remove high frequency noise
    # Retriving Temperature Data
    conn.commit() # Making sure all changes made to database during cycling are commited
    Temp = c.execute('Select MeasT from ' + dataname)   # Getting data from db
    Temp = [item[0] for item in Temp.fetchall()]        # Removing touples

    # Filtering
    b, a = butter(1,w, 'low')                           # Designing Filter
    Tempfill = lfilter(b,a,Temp)                        # Filtering data using low pass filter
    
    # Adding filtered data to db
    for i in range(0,len(Tempfill)):
        c.execute('UPDATE '+dataname+' SET Temp = ? WHERE MeasT = ?',(Tempfill[i]+0.45,Temp[i]))   # Updating last collumn
    conn.commit()                                                                   # Commiting changes
    
    conn.close()            # Closing database connection
    print('Run Complete')   # Confirming that run is finished to console

    return


#---------------------------------------------------------------------------

# Graphical Interface

# Extracting values from database
conn2 = sqlite3.connect('CyclingData.db')    # Connecting to database
c2 = conn2.cursor()                           # Creating cursor for acces

# Extracting standard values from database
profile_time = c2.execute('SELECT Time FROM tempprofile')
profile_time = [int(item[0]) for item in profile_time.fetchall()]
profile_temp = c2.execute('SELECT Temp FROM tempprofile')
profile_temp = [int(item[0]) for item in profile_temp.fetchall()]

tDE = profile_time[0] # Time for denaturation
Tde = profile_temp[0] # Temperature for denaturation step       
tAN = profile_time[1] # Time for annealing step
Tan = profile_temp[1] # Temperature for annealing step
tEX = profile_time[2] # Time for extension step
Tex = profile_temp[2] # Temperature for extesnion step
Nc  = profile_time[3] # Number of cycles
tHS = profile_time[4] # Time for hot start
Ths = profile_temp[4] # Temperatrue for hot start
conn2.close()

root = tk.Tk() # Initialising graphical interface

# Header
logoSogang = tk.PhotoImage(file='Sogang.gif')
h1 = tk.Label(root,image=logoSogang).grid(row=1, column=3)
h2 = tk.Label(root,justify='left',padx=10,text="dPCR Thermalcycling Software V1",font="Verdan 24 bold").grid(row=1, column=1)
h3 = tk.Label(root,justify='left',padx=10,text="by C.D. Ahrberg                                                 2018 ",font="Verdan 16 italic").grid(row=2, column=1)
h3 = tk.Label(root,justify='left',padx=10,text="     ",font="Verdan 24 italic").grid(row=3, column=1)

# Textboxes for entering values
# Hot start
e1 = tk.Entry(root) # Time for hot start in seconds
e1.insert(0, tHS)
e12= tk.Entry(root) # Temperature for Hot start
e12.insert(0, Ths)
# Number of cycles
e2 = tk.Entry(root) # Number of cycles
e2.insert(0, Nc)
# Denaturation
e3 = tk.Entry(root) # Time for Denaturation step in seconds
e3.insert(0, tDE)
e32= tk.Entry(root) # Temperature for Denaturation start
e32.insert(0, Tde)
# Annealing
e4 = tk.Entry(root) # Time for Annealing step in seconds
e4.insert(0, tAN)
e42= tk.Entry(root) # Temperature for Annealing start
e42.insert(0, Tan)
# Extension
e5 = tk.Entry(root) # Time for Extension step in seconds
e5.insert(0, tEX)
e52= tk.Entry(root) # Temperature for Extension start
e52.insert(0, Tex)

# Palcing Textboxes
# Hot start
tk.Label(root,text="Time for hot Start in seconds").grid(row=4, column=0)
e1.grid(row=4, column=1)
tk.Label(root,text="Temperature for hot Start in *C").grid(row=4, column =2)
e12.grid(row=4, column=3)
# Cycles
tk.Label(root,text="Number of cycles").grid(row=5)
e2.grid(row=5, column=1)
# Denaturation
tk.Label(root,text="Time for Denaturation Step in seconds").grid(row=6, column=0)
e3.grid(row=6, column=1)
tk.Label(root,text="Temperature for Denaturation Step in *C").grid(row=6, column =2)
e32.grid(row=6, column=3)
# Annealing
tk.Label(root,text="Time for Annealing Step in seconds").grid(row=7, column=0)
e4.grid(row=7, column=1)
tk.Label(root,text="Temperature for Annealing Step in *C").grid(row=7, column =2)
e42.grid(row=7, column=3)
# Extension
tk.Label(root,text="Time for Extension Step in seconds").grid(row=8, column=0)
e5.grid(row=8, column=1)
tk.Label(root,text="Temperature for Extension Step in *C").grid(row=8, column =2)
e52.grid(row=8, column=3)

tk.Label(root, text="  ").grid(row=9)

# Start Button
button = tk.Button(root, text='Start', width=15, bg='green', command = lambda: writeTprof(e2,e1,e12,e3,e32,e4,e42,e5,e52))
button.grid(row=10, column=1)

tk.Label(root, text="  ").grid(row=11)

# Displaying interface
root.mainloop()
