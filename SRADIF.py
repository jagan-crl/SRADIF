def OTLEMS_NIST(Tsol):
    #Based on given NIST radiation absorption coefficients
    if Tsol > 2900.0e0: Tsol = 2900.0e0
    if Tsol < 300.0e0: Tsol = 300.0e0
    # H2O
    EMSW = -2.3093E-01 - 1.1239E+03*(1/Tsol) \
        + 9.4153E+06*(1/Tsol)**2.0 - 2.9988E+09*(1/Tsol)**3.0 \
        + 5.1382E+11*(1/Tsol)**4.0 - 1.8684E+10*(1/Tsol)**5.0   
    # CO2
    EMSD = 1.8741E+01 - 1.2131E+05*(1/Tsol) \
        + 2.735E+08*(1/Tsol)**2.0 - 1.9405E+11*(1/Tsol)**3.0 \
        + 5.631E+13*(1/Tsol)**4.0 - 5.8169E+15*(1/Tsol)**5.0       
    # CO
    if (Tsol <= 750.0e0):
        EMSMO = 4.7869E+00 - 6.953E-02*Tsol \
        + 2.95775E-04*Tsol**2.0 - 4.25732E-07*Tsol**3.0 \
        + 2.02894E-10*Tsol**4.0
    else:
        EMSMO = 1.0090E+01 - 1.183E-02*Tsol \
        + 4.77530E-06*Tsol**2.0 - 5.87209E-10*Tsol**3.0 \
        - 2.53340E-14*Tsol**4.0    
    # CH4
    EMSME = 6.6334E+00 - 3.5686E-03*Tsol \
        + 1.6682E-08*Tsol**2.0 + 2.5611E-10*Tsol**3.0 \
        - 2.6558E-14*Tsol**4.0
    # HF
    if (Tsol <= 700.0e0):
        EMSHF = 435.74 - 2.7099*Tsol \
        + 6.6504E-3*Tsol**2.0 - 7.4327E-6*Tsol**3.0 \
        + 3.1531E-9*Tsol**4.0
    else:
        EMSHF = 21.322 - 4.3028E-2*Tsol \
        + 3.7464E-5*Tsol**2.0 - 1.5053E-8*Tsol**3.0 \
        + 2.2591E-12*Tsol**4.0
    # CF2O
    if (Tsol <= 700.0e0):
        EMSCF = -3.22504E+01 - 6.02322E-02*Tsol \
        + 1.39465E-03*Tsol**2.0 - 2.79470E-06*Tsol**3.0 \
        + 1.59079E-09*Tsol**4.0
    else:
        EMSCF = 1.05084E+02 - 1.55847E-01*Tsol \
        + 8.96794E-05*Tsol**2.0 - 2.28595E-08*Tsol**3.0 \
        + 2.10067E-12*Tsol**4.0
    return EMSW,EMSD,EMSMO,EMSME,EMSHF,EMSCF

def polyfitter(t,R,numfit,order):
    import numpy as np
    import numpy.polynomial.polynomial as poly
    m = int((numfit-1)/2)
    len_exp_size = np.shape(t)
    len_exp = len_exp_size[0]
    rows = len_exp-2*m
    Rf = np.zeros((rows,1))
    dRfdt = np.zeros((rows,1))   
    for i in range(m,len_exp-m):
        f, stats = poly.polyfit(t[i-m:i+m],R[i-m:i+m],order,full=True)
        Rf[i-m] = poly.polyval(t[i],f)
        fprime = poly.polyder(f)
        dRfdt[i-m]= poly.polyval(t[i],fprime)    
    t_new = t[m:len_exp-m] #s
    t_new = np.transpose(t_new) #ms 
    K = (2/Rf)*dRfdt;
    return Rf,dRfdt,t_new,K

def sphinter(xArray,yArray,xValue):                                                 # this function performs the linear interpolation with experimental data to determine a time scale
    import numpy as np
    for i in range(0,np.size(xArray)-1):
        if (xArray[i] <= xValue) and (xArray[i+1] > xValue):
            ratio = (xValue - xArray[i]) / (xArray[i+1] - xArray[i])
            yValue = yArray[i] + ratio * (yArray[i+1] - yArray[i])
    return yValue

def volToRad(Vol_shell,nlay):                                                       # this function takes the array of shell volumes and generates an array of shell outer radii
    import numpy as np
    R_shell = np.zeros(nlay)
    vBurn = 0.0e0
    for jlay in range(0,nlay):
        vBurn = vBurn + Vol_shell[jlay] 
        R_shell[jlay] = ((vBurn * 0.75e0 / np.pi)**(1.0e0/3.0e0))
    return R_shell
    
def spherical(reactants,scheme,pres,temp,nlay,volu,atol,\
              gtol,ttol,poly_switch,ign):
    import sys
    import cantera as ct
    import numpy as np
    import numpy.polynomial.polynomial as poly
    
    gas = ct.Solution(scheme)
    pres=pres*ct.one_atm
    
    
    #################################################
    #     Initialize RSPHWRK
    #################################################
    gas.TPX = temp, pres, reactants                                                 # set initial gas state
    
    ##################################
    # Initialization of the spherical shells with equal width but varying volume
    ##################################
    if ign == 0:                                                                    # all shells same width
        dR1 = (volu * 0.75e0 / np.pi)**(1.0e0/3.0e0) / (float(nlay)-1)
    else:                                                                           # first shell varying width (compensate for initial flame kernel size)
        vol = volu - ign
        dR1 = (vol * 0.75e0 / np.pi)**(1.0e0/3.0e0) / (float(nlay)-1)               
        dR2 = (ign * 0.75e0 / np.pi)**(1.0e0/3.0e0)
    
    ##################################
    #  Arrays which represent state of each shell at the current analysis time
    ##################################
    T_shell = np.array([temp] * nlay)                                               # gas temperature
    R_shell = np.array([dR1] * nlay)
    Urad_shell = np.array([0.0] * nlay)
    if ign > 0:
        R_shell[0] = dR2
    for i in range(0,nlay):
        if i >= 1:
            R_shell[i] = R_shell[i-1]+dR1
    Vol_shell = np.array([4.0e0/3.0e0 * np.pi * R_shell[0]**3.0e0])                 # [cm^3] array of shell volumes
    for i in range(1,nlay):
        V1 = 4.0e0/3.0e0 * np.pi * R_shell[i]**3.0e0
        V2 = 4.0e0/3.0e0 * np.pi * R_shell[i-1]**3.0e0
        Vol_shell = np.append(Vol_shell, V1-V2)
    Gamma_shell = np.array([gas.cp/gas.cv]*nlay)                                    # unburned gas specific heat ratio cp/cv   
    Species_shell = np.array([gas.X] * nlay)                                        # array of shell species
    Species_shell_Y = np.array([gas.Y] * nlay)                                        # array of shell species

    ##################################      
    # Arrays for tracking values from before radiative cooling is applied
    ##################################
    T_old = np.zeros(nlay)
    Vol_old = np.zeros(nlay)
    Gamma_old = np.zeros(nlay)
    R_old = [0.0]
    
    ##################################
    # Arrays for saving temperature profile diagnostics
    ##################################
    T_test = [0.0]
    R_test = [0.0]
    xCO_test = [0.0]
    xCO2_test = [0.0]
    xHF_test = [0.0]
    yCO_test = [0.0]
    yCO2_test = [0.0]
    yHF_test = [0.0]
    dt_test = [0.0]
    iter_test = [0.0]
    check2 = 0
    
    ##################################
    # Arrays which represent the history of a particular variable over time
    ##################################
    delt_his = [0.0]
    P_his = [pres]                                                                  # should remain constant pressure
    T_his = [temp]
    T_core = [temp]
    R_his = [0.0]
    dR = [0.0]
    ub = [0.0]
    time_his = [0.0]                                                                # [s] time history
    
    
    ##################################
    #  Initialize experiment data array
    ################################## 
    
    # Experimental data input
    expdata = np.loadtxt('globalOutput.dat', encoding="utf8")                       # Specify filename for importing experimental data
    
    ##################################
    # Unit conversions
    ##################################

    #If using SLTORC globalOutput File
    texp = expdata[:,0] * 1000                                                      # [s] to [ms]
    rexp = expdata[:,3] * 100                                                       # [m] to [cm]

    ##################################
    #  Interpolating polynomial for input data
    ##################################
    
    if poly_switch == 1:
        cutoff = 0                                                                  # how many initial data points are neglected from each side
        polyorder = 4                                                               # order of interpolating polynomial                                          

        texp_fit = texp[cutoff:]
        rexp_fit = rexp[cutoff:]
        coefs = poly.polyfit(texp_fit, rexp_fit, polyorder)
        
        texp = np.linspace(texp_fit[0], texp_fit[-1], num=3800)
        rexp = poly.polyval(texp, coefs)
   
    ##################################
    #  Read some indices for radiation calculations
    ##################################
    
    kH2O = gas.species_index('H2O')
    kCO2 = gas.species_index('CO2')
    kCO = gas.species_index('CO')
    kCH4 = gas.species_index('CH4')
    kHF = gas.species_index('HF')
    kCF2O = gas.species_index('CF2O')
    Sigma = 5.6704e-8                                                               # [(J/s)/(m^2*K^4)]
    #     In cgs units the Stefan-Boltzmann constant: erg*cm^-2*s^-1*K-4
    Rideal = 82.05736e-3
    #     ideal gas constant = 82.05736 cm^3*atm/(K*mol) 
    #     ideal gas constant = 0.08205736 m^3*atm/(K*kmol)
        
    ##################################
    #  Initiate printing profiles at specified timesteps
    ##################################
    
    nDeltat = 1.0
    writeDeltat = 1.0                                                               # Time interval for writing data to profiles file [ms]
    writeNow = 0.0
    
    profiles=open('profiles.dat','w')                                               # Specify name of profiles output file
    profiles.write("#      r(cm)        T(K) u_rad(cm/s)        xH2O        xCO2         xCO        xCH4         xHF       xCF2O\n")
    profiles.flush()
    
    #################################################
    #  Start Computations
    #################################################
    print('\n---- INITIAL GAS STATE')
    gas()                                                                           # prints intitial gas properties to screen
    for ilay in range(0, nlay):
        print(ilay)                
        #################################################
        #   Dissociation calculation for inner shells
        #################################################
        for jlay in range(0,ilay):
            gas.TPX = T_shell[jlay], P_his[ilay], Species_shell[jlay]               # create gas object with specifications of burned gas shell
            vini = gas.v                                                            # save specific volume before and after to determine volume change
            gas.equilibrate('HP')
            vf = gas.v;
            T_shell[jlay] = gas.T                                                   # burned gas temperature
            Gamma_shell[jlay] = gas.cp/gas.cv                                       # layer specific heat ratio cp/cv
            Vol_shell[jlay] = (Vol_shell[jlay])*(vf/vini)                           # determine change in volume by the ratio of specific volumes
            Species_shell[jlay] = gas.X                                             # burned species composition
            Species_shell_Y[jlay] = gas.Y                                           # burned species composition   
        
        #################################################
        #  Combustion equilibrium of the burning layer
        #################################################
        gas.TPX = T_shell[ilay], P_his[ilay], Species_shell[ilay]                   # create gas object with specifications of burned gas shell
        vini = gas.v                                                                # save specific volume before and after to determine volume change
        gas.equilibrate('HP')
        vf = gas.v
        T_shell[ilay] = gas.T                                                       # burned gas temperature
        Gamma_shell[ilay] = gas.cp/gas.cv                                           # layer specific heat ratio cp/cv
        UVOL = Vol_shell[ilay]                                                      # used for initial delta-time guess
        Vol_shell[ilay] = (Vol_shell[ilay])*(vf/vini)                               # determine change in volume by the ratio of specific volumes
        Species_shell[ilay] = gas.X                                                 # burned species composition
        Species_shell_Y[ilay] = gas.Y                                               # burned species composition

        delt = 0.0e0
        
        Lrad = True
        if Lrad==1:
        #################################################
        #   Radiation Calculation
        #################################################
        
            ##################################
            #  Pick an initial delt guess:
            ##################################
            Su = (P_his[ilay]/ct.one_atm) * 5.0e0  #cm/sec
            if (ilay > 0) and (Su > 0):
                delt = UVOL/(4.0e0*np.pi*R_his[ilay]**2.0e0*Su)
            else:
                delt = 0.0e0
             
            ##################################
            #  If experiment data used, then save initial state for time iteration
            ##################################
            for jlay in range (0,nlay):                                             # save pre-radiation cooling values
                Vol_old[jlay] = Vol_shell[jlay]
                T_old[jlay] = T_shell[jlay]
                Gamma_old[jlay] = Gamma_shell[jlay]

            ##################################
            # OTL calculation of RADENERGY  #K/sec
            ##################################
            Qrad = np.zeros(ilay+1)
            dTdt_rad = np.zeros(ilay+1)
            Crad = 1.0                                                              # radiation scaling coefficient
            for jlay in range(0,ilay+1):
                Tsol = T_shell[jlay]
                (EMSW,EMSD,EMSMO,EMSME,EMSHF,EMSCF) = OTLEMS_NIST(Tsol)             # use polynimials
                SumE = 0.0e0
                
                SumE = SumE + EMSW * Species_shell[jlay][kH2O] \
                    + EMSD * Species_shell[jlay][kCO2] \
                    + EMSMO * Species_shell[jlay][kCO] \
                    + EMSME * Species_shell[jlay][kCH4] \
                    + EMSHF * Species_shell[jlay][kHF] \
                    + EMSCF * Species_shell[jlay][kCF2O]                            # add values extracted from coefficients, by species amount
                
                EMIST = SumE * (P_his[ilay]/ct.one_atm)                             # [1/m] pressure should be in atm
                
                gas.TPX = T_shell[jlay], P_his[ilay], Species_shell[jlay]           # prepare CP value in units that will cancel with Qdot to extract temperature change per time change
                CP = gas.cp_mole                                                    # Cp in molar basis [J/kmol/K]
                CP = CP * (P_his[ilay]/ct.one_atm)/(Rideal*T_shell[jlay])           # Cp*(Kmol/m^3) = [J/m^3/K]
                
                RADENERGY = 4.0e0 * EMIST * Sigma * Crad \
                    * (T_shell[jlay]**4.0e0 - temp**4.0e0)                          # [(J/s)/m^3]
                Qrad[jlay] = RADENERGY                                              # [W/m^3] Qdot must be considered separately for the radiation case
                RADENERGY = RADENERGY / CP
                dTdt_rad[jlay] = RADENERGY                                          # [K/s]

            
            ##################################
            #  Cooling time iteration for radiation
            ##################################
            Vtemp = 0.0e0
            for jlay in range(0,ilay+1):
                Vtemp = Vtemp + Vol_shell[jlay]
            R_old.append((Vtemp * 0.75e0 / np.pi)**(1.0e0/3.0e0))                   # find current burned gas radius
            
            R_shell_old = volToRad(Vol_shell,nlay)
            
            Iter = 0
            IterMax = 1000
            rlx = 0.7                                                               # delta scaling coefficient
            deltOld = delt + ttol * 10.0e0                                          # initial delta-time old guess for checking tolerances
            while (abs(delt-deltOld) >= ttol):                                      # delta-time determination loop
                Iter = Iter + 1
                if Iter > 1:
                    delt = deltOld * rlx + delt * (1.0 - rlx)
                
                ##################################
                #  Cooling of each layer
                ##################################

                for jlay in range(0,ilay+1):                                        # find temperature of each shell after radiation effects are applied
                    Trad = T_old[jlay] - (dTdt_rad[jlay])*delt
                    if Trad <= 0:
                        print('oof')
                    Vol_shell[jlay] = Vol_old[jlay] * Trad / T_old[jlay]            # move dissociation here, v/v
                    T_shell[jlay] = Trad 
                    gas.TPX = T_shell[jlay], P_his[ilay], Species_shell[jlay]
                    Gamma_shell[jlay] = gas.cp/gas.cv                               # layer specific heat ratio cp/cv

                #################################################
                #   Iteratively determine delt
                #################################################
                # linear interpolation performed with radius
                # Hold pressure constant (no contact with chamber walls)          
                Vtemp = 0.0e0
                for jlay in range(0,ilay+1):
                    Vtemp = Vtemp + Vol_shell[jlay]
                Rtemp = (Vtemp * 0.75e0 / np.pi)**(1.0e0/3.0e0)                     # temporary variable for radius after radiation-induced volume decrease
                
                deltOld = delt
                if (R_his[ilay] >= rexp[0]) and (Rtemp <= rexp[rexp.size-1]):
                    if (Rtemp <= R_his[ilay]):                                      # if volume increasing with radiation
                        if (deltOld <= ttol):
                            out.write("Error in Spherical: radius not decreasing.")
                            out.write(R_his[ilay],Rtemp)
                            out.close() 
                            exit()
                        delt = (deltOld - ttol) / 2.0
                    else:                                                           # linearly interpolate with experimental data to determine delta-time
                        if poly_switch == 1:
                            time1 = sphinter(r_new,texp,R_his[ilay])
                            time2 = sphinter(r_new,texp,Rtemp)
                        else:
                            time1 = sphinter(rexp,texp,R_his[ilay])
                            time2 = sphinter(rexp,texp,Rtemp)                            
                        delt = (time2 - time1) / 1000.0e0                           # convert [ms] to [s]
                        if (delt <= 0.0e0):
                            out.write("Error in exppt.txt: radius not decreasing & delt < 0.")
                            out.write(R_his[ilay],time1,Rtemp,time2)
                            out.close() 
                            exit()
                if (Iter >= IterMax):
                    out.write("Cooling time iteration exceed ",IterMax," times.")
                    out.close()
                    sys.exit()
                    
            delt_his.append(delt)                                                   # save delta-time

        
        #################################################
        #   Calculates flame radius
        #################################################        
        P_his.append(P_his[ilay])                                                   # pressure assumed constant throughout this process
        
        vBurn = 0.0e0
        for jlay in range(0,ilay+1):
            vBurn = vBurn + Vol_shell[jlay]      
        R_his.append((vBurn * 0.75e0 / np.pi)**(1.0e0/3.0e0))                       # find current burned gas sphere radius
        
        T_his.append(T_shell[nlay-1])                                               # temperature of unburned gas
        T_core.append(T_shell[0])                                                   # temperature of the first shell
        
        R_shell = volToRad(Vol_shell,nlay)                                          # update array of all current layer radii
                
        #################################################
        #   Calculates inflow velocity
        #################################################
        
        if ilay==0:
            ub[ilay] = 0
        else:  
            ub.append((R_old[ilay+1]-R_his[ilay+1])/(delt))                         # delta-time
            dR.append(R_old[ilay+1]-R_his[ilay+1])                                  # delta-radius
            for jlay in range(0,nlay):
                Urad_shell[jlay] = (R_shell[jlay]-R_shell_old[jlay])/delt
        
        if (R_his[ilay] >= rexp[0]) and (Rtemp <= rexp[rexp.size-1]):               # append current simulation time
            time_his.append(time2/1000)
        else:
            time_his.append(0)
                
        #################################################
        #   Writing profiles to output file
        ################################################# 

        tnow = 1000*time_his[ilay+1] #ms
        
        if tnow >= writeNow:
            profiles.write("#%1.0f - t=%8.6fms\n"%(nDeltat,tnow))
            for n in range(0,nlay):
                profiles.write("%12.4e"%(R_shell[n]))
                profiles.write("%12.4e"%(T_shell[n]))
                profiles.write("%12.4e"%(Urad_shell[n]))
                profiles.write("%12.4e"%(Species_shell[n][kH2O]))
                profiles.write("%12.4e"%(Species_shell[n][kCO2]))
                profiles.write("%12.4e"%(Species_shell[n][kCO]))
                profiles.write("%12.4e"%(Species_shell[n][kCH4]))
                profiles.write("%12.4e"%(Species_shell[n][kHF]))
                profiles.write("%12.4e\n"%(Species_shell[n][kCF2O]))
            profiles.write("\n\n")
            profiles.flush()
            nDeltat = nDeltat + 1
            writeNow = writeNow + writeDeltat
  
        
    #################################################
    # Process output data
    #################################################
    
    # Initialize polynomial fit
    numfit = nlay/10+1
    order = 2
    
    Rf_sol,dRfdt_sol,t_new,K = polyfitter(time_his,R_his,numfit,order)
    
    m = int((numfit-1)/2)
    len_exp_size = np.shape(time_his)
    len_exp = len_exp_size[0]
    rows = len_exp-2*m
    ub_scale = np.zeros((rows,1))
    for i in range(m,len_exp-m):
        ub_scale[i-m] = ub[i]

    #################################################
    # Export data
    #################################################
    
    out = np.column_stack((1000*t_new,Rf_sol,K,dRfdt_sol,-ub_scale,dRfdt_sol+ub_scale))
    np.savetxt("Output.txt",out)                                                    # Specify name of output data file
    
