def calc_EIS(Ts, SLP, T700, z700):
    '''
    Input:
           Ts: surface temperature (K)
           SLP: sea level pressure (Pa)
           T700: temperature at 700 hPa (K)
           z700: height at 700 hPa (m)
    Output: 
           EIS: estimated inversion strengh (K)
           LTS: lower troposheric strength (K)

    # Notes:
    #  - surface relative humidity is assumed to be constant
    #  - potential temperature is calcualted with respect to the 1000 hPa
    #    reference level
    #  - traditionally, surface air temperature is used in place of surface
    #    temperature.
    #  - physical constants are from Curry and Webster, Appendix B
    #
    # Physical Constants:
    #  - Gas constant for water vapor = 461 J/K/kg
    #  - Latent heat of vaporization at 273.15 K = 2.5*10^6 J/kg
    #  - Gas constant for dry air = 287 J/K/kg
    #  - Specific heat at constant pressure = 1004 J/K/kg
    #
    #References:
    #
    #Curry, J. A., & Webster, P. J. (1998). Thermodynamics of Atmospheres and 
    #Oceans. Elsevier.
    #
    #Georgakakos, K. P., & Bras, R. L. (1984). A hydrologically useful station
    #precipitation model: 1. Formulation. Water Resources Research, 20(11), 
    #1585-1596, https://doi.org/10.1029/WR020i011p01585
    #
    #Wood, R., & Bretherton, C. S. (2006). On the relationship between 
    #stratiform low cloud cover and lower-tropospheric stability. Journal of 
    #Climate, 19(24), 6425-6432, https://doi.org/10.1175/JCLI3988.1

    '''

    #calculate the LCL height
    es=6.11*np.exp(((2.5*10**6)/461)*((1/273.15)-(1/Ts))) #Clausius-Clapeyron Equation (Equation 4.23 from Curry and Webster)
    e=es*80/100 #assume RH is 80%, as in WB06
    Td=((1/273.15)-(461/(2.5*10**6))*np.log(e/6.11))**(-1) #form of Clausius-Clapeyron Equation
    p_LCL=SLP*(((Ts-Td)/223.15)+1)**(-3.5) #Equation (12) from Georgakakos and Bras (1984)
    T_LCL=Ts*(((Ts-Td)/223.15)+1)**(-1) #Equation (13) from Georgakakos and Bras (1984)
    LCL=((287*(Ts+T_LCL)/2)/9.8)*np.log(SLP/p_LCL) #Hypsometric equation

    #calculate LTS
    theta700=T700*(1000/700)**(287/1004)
    theta_slp=Ts*(1000./(SLP*0.01))**(287/1004)
    LTS=theta700-theta_slp

    #calculate the moist adaiabtic lapse rate at 850 hPa
    T_bar=(Ts+T700)/2 #approximation of the temperature at 850 hPa, following WB06
    es_bar=6.11*np.exp(((2.5*10**6)/461)*((1/273.15)-(1/T_bar))) #Clausius-Clapeyron Equation (Equation 4.23 from Curry and Webster)
    qs_bar=0.622*es_bar/(850-es_bar) #Equation 4.37 of Curry and Webster
    gamma_m=(9.8/1004)*(1-((1+(((2.5*10**6)*qs_bar)/(287*T_bar)))/(1+((((2.5*10**6)**2)*qs_bar)/(1004*461*(T_bar**2)))))) #Equation (5) of WB06
    
    EIS=LTS-gamma_m*(z700-LCL); #Equation (4) of WB06

    return EIS, LTS


