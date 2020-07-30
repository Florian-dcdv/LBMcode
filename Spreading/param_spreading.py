testWe=[10,30,50,100]
testD=[40,50,60]
gmin=1000
gmax=0
for We in testWe:
    for D in testD:
        a=0.50
        Tc=a/(2.65018*4)
        #We=52
        tauv=0.6
        nu=(tauv-0.5)/3
        drho=0.431-0.00149     #density difference
        #drho=0.27582-0.034
        #gamma = 1.265e-2*a + 3.892e-3        #surface tension
        gamma=8.871e-3
        H_D=2
        u=(gamma*We/(drho*D))**0.5
        Re=u*D/nu
        Oh=We**0.5/Re
        H=H_D*D
        g=u**2/(2*(H-D))    #gravity
        Fb=g*drho           #body force
        print("a=",a,"| Tc=",Tc)
        print()
        print("Dimensionless numbers:")
        print("We =",We)
        print("Re =",Re)
        print("Oh =",Oh)
        print("H/D =",H_D)
        print()
        print("Fixed parameters:")
        print("D =",D)
        print("rho =",drho)
        print("gamma =",gamma)
        print("tauv=",tauv)
        print("viscosity=",nu)
        print()
        print("Calculated values:")
        print("H =",H)
        print("u =",u)
        print("g=",g)
        print("Fb=",Fb)
        print()
        print("#############################")
        print()
        if g<gmin:
            gmin=g
        elif g>gmax:
            gmax=g
        print("gmin=",gmin,"gmax=",gmax)


