! FORTRAN MODULE FOR CALCULATING THE CLUSTER FORMATION EFFICIENCY (GLOBAL VERSION)
! Copyright (C) 2012  Diederik Kruijssen
!
! NAME:
!       CFEmod (module)
!       F_CFE (main function)
!
! PURPOSE:
!       Calculate the fraction of star formation occurring in bound clusters,
!       i.e. the cluster formation efficiency or CFE
!
! TERMS OF USE:
!       If you use this routine while preparing a paper, please cite the paper
!       in which this model was presented:
!       Kruijssen, J. M. D., 2012, MNRAS 426, 3008
!
! CALLING SEQUENCE:
!       f_cfe(surfg,qT,omega)
!
! INPUT PARAMETERS:
!       NOTE: ALL INPUT SHOULD BE IN SI UNITS!
!       surfg - gas surface density
!       qT - Toomre (1964) Q parameter (if unknown set to 1.5)
!       omega - angular velocity (if unknown use function: omega=f_omega(surfg))
!
! OPTIONAL INPUT IS SET IN THE FILE "PARAMETERS"
!       sflaw - star formation law: 0=Elmegreen (2002, default)
!                                   1=Krumholz & McKee (2005)
!       qvir - giant molecular cloud virial parameter (default=1.3)
!       tsn - time of the first supernova (default=3.*Myr)
!       tview - time at which CFE is determined (default=10.*Myr)
!       surfGMC - giant molecular cloud surface density (default=100.*Msun/pc**2.)
!       ecore - maximum (protostellar core) star formation efficiency (default=.5)
!       beta0 - turbulent-to-magnetic pressure ratio (default=1.e10)
!       radfb - feedback: 0=supernova, 1=radiative, 2=both (default=0)
!
! OUTPUT:
!       cfearray - array with four elements, the 1st being the CFE, the 2nd the
!       naturally bound fraction of star formation, the third the fraction of
!       bound star formation that survives the cruel cradle effect, and the 4th
!       the fraction of all star formation that survives the cruel cradle effect
!
! EXAMPLE IS GIVEN IN THE FILE testCFE.F90:
!       Calculate the CFE in the solar neighbourhood
!
!       PROGRAM testCFE
!           use CFEmod
!           REAL pc,msun,myr,surfMW,qMW,omegaMW,cfe(1:4)
!           pc=3.1e16
!           Msun=2.e30
!           Myr=3.16e13
!           surfMW=9.3*msun/pc**2.
!           qMW=2.0
!           omegaMW=0.026/Myr
!           cfe=f_cfe(surfMW,qMW,omegaMW)
!           PRINT*,cfe(1)*100. 
!       END PROGRAM
!
!       > g95 f_cfe.f90 testCFE.f90 -o testCFE
!       > ./testCFE
!          6.702944
!
! MODULE STRUCTURE: 
!       Auxiliary functions come first, the actual CFE function is located at
!       the bottom of this file
!
! REVISION HISTORY:
!       Written by Diederik Kruijssen, August 2012
!

MODULE CFEmod
    CONTAINS

FUNCTION f_tff(rho) result(tff) !free-fall time as a function of density
    REAL, INTENT(in) :: rho
    REAL :: G,pi,tff
    G=6.67e-11 !gravitational constant
    pi=3.14159265358979323846 !pi
    tff=sqrt(3.*pi/32./G/rho) !free-fall time
END FUNCTION

FUNCTION f_fgmc(surfg) result(fgmc) !molecular gas fraction as a function of surface density, from Krumholz & McKee (2005)
    REAL, INTENT(in) :: surfg
    REAL :: msun,pc,surfg2,fgmc
    msun=1.989e30 !solar mass in kg
    pc=3.086e16 !parsec in meters
    surfg2=surfg/100./msun*pc**2. !surface density in units of 100 Msun/pc**2
    fgmc=1./(1.+0.025*surfg2**(-2.)) !molecular gas fraction
END FUNCTION

FUNCTION f_phi(surfg) result(phi) !ratio between cloud pressure and mid-plane pressure as a function of surface density, from Krumholz & McKee (2005)
    REAL, INTENT(in) :: surfg
    REAL :: fgmc,phi
    fgmc=f_fgmc(surfg) !molecular gas fraction
    phi=(10.-8.*fgmc) !ratio between cloud pressure and mid-plane pressure
END FUNCTION

FUNCTION f_mach(surfg,qT,omega) result(mach) !Mach number as a function of surface density, angular velocity and Toomre Q, from Krumholz & McKee (2005)
    REAL, INTENT(in) :: surfg,qT,omega
    REAL :: msun,pc,myr,omega0,phi,mach
    msun=1.989e30 !solar mass in kg
    pc=3.086e16 !parsec in meters
    myr=1.e6*86400.*365.25 !million years in seconds
    omega0=omega*myr !angular velocity in units of 1/Myr
    phi=f_phi(surfg) !ratio between cloud pressure and mid-plane pressure, see above
    mach=2.82*phi**.125*qT/omega0*surfg/100./msun*pc**2. !Mach number
END FUNCTION

FUNCTION f_sigrho(mach,beta0) result(sig) !dispersion of overdensity PDF as a function of Mach number and magnetic pressure ratio, based on e.g. Padoan & Nordlund (2011)
    REAL, INTENT(in) :: mach,beta0
    REAL :: b,sig
    b=0.5 !constant
    sig=SQRT(LOG(1.+3.*b**2.*mach**2.*beta0/(beta0+1.))) !dispersion
END FUNCTION

FUNCTION f_xcrit(qvir,mach) result(x) !critical overdensity for SF in the KM05 sSFR_ff as a function of Mach number and GMC virial ratio, from Krumholz & McKee (2005)
    REAL, INTENT(in) :: qvir,mach
    REAL :: pi,phix,x
    pi=3.14159265358979323846 !pi
    phix=1.12 !constant
    x=pi**2.*phix**2./15.*qvir*mach**2. !critical overdensity
END FUNCTION

FUNCTION f_sfrff(qvir,mach,beta0,ecore,sflaw) result(f) !specific star formation rate per free-fall time (sSFR_ff) for Elmegreen (2002, sflaw=0) or Krumholz & McKee (2005, sflaw=1)
    REAL, INTENT(in) :: qvir,mach,beta0,ecore
    REAL :: phit,xcrit,sigrho,f
    INTEGER, INTENT(in) :: sflaw
    phit=1.91 !constant
    xcrit=f_xcrit(qvir,mach) !critical overdensity for star formation in the KM05 model, see above
    sigrho=f_sigrho(mach,beta0) !overdensity PDF dispersion, see above
    IF(sflaw.EQ.1) THEN
        f=.5*ecore/phit*(1.+erf((-2.*alog(xcrit)+sigrho**2.)/(2.**1.5*sigrho))) !specific star formation rate per free-fall time
    ELSE
        f=0.012 !specific star formation rate per free-fall time
    ENDIF
END FUNCTION

FUNCTION f_dpdx(x,mulnx,sig) result(dpdx) !overdensity PDF of the ISM as a function of overdensity x and its logarithmic mean and dispersion
    REAL, INTENT(in) :: x,mulnx,sig
    REAL :: pi,dpdx
    pi=3.14159265358979323846 !pi
    dpdx=1./(sqrt(2.*pi*sig**2.)*x)*exp(-.5*(alog(x)-mulnx)**2./sig**2.) !overdensity PDF
END FUNCTION

FUNCTION f_rho0(qT,omega) result(rho0) !mid-plane ISM density as a function of Toomre Q and angular velocity, assuming hydrostatic equilibrium (cf. Krumholz & McKee 2005)
    REAL :: qT,omega
    REAL :: G,pi,phiP,rho0
    G=6.67e-11 !gravitational constant
    pi=3.14159265358979323846 !pi
    phiP=3. !constant
    rho0=phiP*omega**2./(pi*qT**2.*G) !mid-plane ISM density
END FUNCTION

FUNCTION f_fstar(surfg,qT,omega,x,ecore,beta0,qvir,tsn,tview,surfGMC,sflaw,radfb) result(fstar) !naturally bound fraction of star formation
    REAL, INTENT(in) :: surfg,qT,omega,x,ecore,beta0,qvir,tsn,tview,surfGMC
    REAL :: G,pi,sigSB,c,phifb,kappa0,psi,phitrap,surffb,mach,sfrff,rhog,tff,efb,einc,efbrad,epsilons(1:4),fstar
    INTEGER, INTENT(in) :: sflaw,radfb

    G=6.67e-11 !gravitational constant
    pi=3.14159265358979323846 !pi
    sigSB=5.67e-8 !Stefan-Boltzmann constant
    c=299792458. !speed of light
    phifb=1.6e-5 !feedback efficiency
    kappa0=2.4e-5 !opacity constant
    psi=.3 !light-to-mass ratio
    phitrap=.2 !trapping ratio

    surffb=MAX(surfGMC,surfg) !surface density on which radiative feedback acts
    mach=f_mach(surfg,qT,omega) !Mach number, see above
    sfrff=f_sfrff(qvir,mach,beta0,ecore,sflaw) !specific star formation rate per free-fall time, see above
    rhog=x*f_rho0(qT,omega) !mid-plane ISM density, see above
    tff=f_tff(rhog) !free-fall time, see above

    !Star Formation Efficiencies
    IF(radfb.EQ.0.OR.radfb.EQ.2) THEN
        efb=0.5*sfrff*tsn/tff*(1.+sqrt(1.+2.*pi**2.*G**2.*tff*qT**2.*surfg**2./(phifb*sfrff*tsn**2.*omega**2.*x))) !SN feedback
    ELSE
        efb=1.
    ENDIF
    einc=sfrff*tview/tff !star formation is incomplete/still ongoing
    IF(radfb.GT.0) THEN
        efbrad=2.*sigSB/(phitrap*kappa0**2.*psi*surffb**3.)*(sqrt(1.+2.*pi*c*G*phitrap*kappa0**2.*surffb**4./(1.*sigSB))-1.) !radiative feedback
    ELSE
        efbrad=1.
    ENDIF
    epsilons=[ecore,efb,efbrad,einc] !SFEs for [maximum,SNfeedback,radiativefeedback,incomplete]
    fstar=MIN(epsilons(1),epsilons(2),epsilons(3),epsilons(4)) !local SFE is the minimum of those
END FUNCTION

FUNCTION f_integrate(xsurv,mulnx,sig,surfg,qT,omega,ecore,beta0,qvir,tsn,tview,surfGMC,cce,sflaw,radfb) result(frac) !obtain fractions from integrating the overdensity PDFs
    REAL, INTENT(in) :: xsurv,mulnx,sig,surfg,qT,omega,ecore,beta0,qvir,tsn,tview,surfGMC
    REAL :: xmin1,xmax,xarr(1:1000),f1,dx,xg,fstar,bound,integral,dpdx,f2,frac
    INTEGER, INTENT(in) :: cce,sflaw,radfb
    INTEGER :: nx,ix

    nx=1000 !number of integration steps (checked to be sufficient for convergence)
    xmin1=EXP(mulnx-5.*sig) !minimum overdensity
    xmax=EXP(mulnx+10.*sig) !maximum overdensity
    IF(cce.GT.0.AND.xsurv.LT.xmin1) xmin1=xsurv !if calculating the cruel cradle effect and critical overdensity below minimum, then adjust
    IF(cce.GT.0.AND.xsurv.GT.xmax) xmax=xsurv !if calculating the cruel cradle effect and critical overdensity below maximum, then adjust
    DO ix=1,nx
        xarr(ix)=xmin1*(xmax/xmin1)**((ix-0.5)/nx) !integration array
    ENDDO
    f1=0. !initialize integral
    DO ix=1,nx !denominator integral
        dx=xarr(ix)*((xmax/xmin1)**(1./(2.*nx))-(xmax/xmin1)**(-1./(2.*nx))) !step size
        xg=xarr(ix) !overdensity
        fstar=f_fstar(surfg,qT,omega,xg,ecore,beta0,qvir,tsn,tview,surfGMC,sflaw,radfb) !local SFE
        bound=fstar/ecore !local bound fraction
        IF(cce.EQ.0) bound=1. !if not calculating the cruel cradle effect but the naturally bound fraction of SF, the denominator should contain all SF
        IF(cce.EQ.2) bound=1. !if calculating the cruel cradle effect with respect to all SF, the denominator should contain all SF
        integral=bound*fstar*xarr(ix) !integral part 1
        dpdx=f_dpdx(xg,mulnx,sig) !overdensity PDF, i.e. integral part 2
        f1=f1+integral*dpdx*dx !integral
    ENDDO
    
    IF(cce.GT.0) THEN !if calculating the cruel cradle effect set minimum overdensity to critical overdensity
        xmin2=xsurv
    ELSE
        xmin2=xmin1
    ENDIF
    DO ix=1,nx
        xarr(ix)=xmin2*(xmax/xmin2)**((ix-0.5)/nx) !integration array
    ENDDO
    f2=0. !initialize integral
    DO ix=1,nx !numerator integral
        dx=xarr(ix)*((xmax/xmin2)**(1./(2.*nx))-(xmax/xmin2)**(-1./(2.*nx))) !step size
        xg=xarr(ix) !overdensity
        fstar=f_fstar(surfg,qT,omega,xg,ecore,beta0,qvir,tsn,tview,surfGMC,sflaw,radfb) !local SFE
        bound=fstar/ecore !local bound fraction
        IF(cce.EQ.2) bound=1.  !if calculating the cruel cradle effect with respect to all SF, the numerator should contain all SF
        integral=bound*fstar*xarr(ix) !integral part 1
        dpdx=f_dpdx(xg,mulnx,sig) !overdensity PDF, i.e. integral part 2
        f2=f2+integral*dpdx*dx !integral
    ENDDO
    
    IF(f1.EQ.0) THEN
        frac=0.
    ELSE
        frac=f2/f1 !numerator divided by denominator
    ENDIF
END FUNCTION

FUNCTION f_phit(qvir,x) result(phit) !ratio of encounter timescale to energy dissipation timescale as a function of cloud virial ratio and overdensity
    REAL, INTENT(in) :: qvir,x
    REAL :: phit
    phit=3.1*SQRT((qvir/1.3)*(x/1.e4)) !ratio of encounter timescale to energy dissipation timescale
END FUNCTION

FUNCTION f_phiad(qvir,x) result(phiad) !adiabatic correction as a function of cloud virial ratio and overdensity
    REAL, INTENT(in) :: qvir,x
    REAL :: phit,phiad
    phit=f_phit(qvir,x) !ratio of encounter timescale to energy dissipation timescale
    phiad=EXP(-2.*phit) !adiabatic correction
END FUNCTION

FUNCTION f_xcce(surfg,qT,omega,surfGMC,qvir,tview) result(xfit) !critical overdensity to remain bound despite the cruel cradle effect (also see Kruijssen et al. 2011)
    REAL, INTENT(in) :: surfg,qT,omega,surfGMC,qvir,tview
    REAL :: pi,eta,g_close,phish,rh2r2av,f,xmin,xmax,xfit,accuracy,xfit0,xarr(1:101),xarr2(1:101),phiad,diff,diffx
    INTEGER :: nx,niter,itermax,ix,ixfit
    
    pi=3.14159265358979323846 !pi
    eta=2.*1.305*3.*pi/64. !for Plummer
    g_close=1.5 !close encounter correction
    phish=2.8 !higher-order energy loss correction
    rh2r2av=.25 !for Plummer
    f=0.7 !fraction of injected energy that is used for unbinding the region
    
    !SOLVE implicit relation for x_cce
    xmin=1.e-4 !minimum x
    xmax=1.e8 !maximum x
    nx=101 !length of x array
    niter=0 !number of elapsed iterations
    itermax=10 !maximum number of iterations
    xfit=xmin**2./xmax !initialisation of fitted x
    accuracy=1.e-6 !desired logarithmic accuracy
    xfit0=xfit*accuracy !initialisation of previously fitted x
    
    DO WHILE(ABS(LOG10(xfit/xfit0)).GT.accuracy.AND.niter.LT.itermax) !while iteration does not give desired convergence, do
        xfit0=xfit !previously fitted x
        DO ix=1,nx !for all x
            xarr(ix)=10.**((ix-1)/(nx-1.)*LOG10(xmax/xmin)+LOG10(xmin)) !x array
            phiad=f_phiad(qvir,xarr(ix)) !adiabatic correction, see above
            xarr2(ix)=124.*g_close*phish*omega*surfGMC*f*phiad*tview/qT/surfg/SQRT(pi) !right-hand side of equation
        ENDDO
        diff=1.e30
        DO ix=1,nx !for all x
            diffx=ABS(xarr(ix)-xarr2(ix))
            IF(diffx.LT.diff.AND.ix.GT.1) THEN
                diff=diffx
                ixfit=ix !index where x equals right-hand side of equation
            ENDIF
        ENDDO
        xfit=xarr(ixfit) !solution for x_cce
        xmin=xarr(ixfit-1) !new minimum x
        xmax=xarr(ixfit+1) !new maximum x
        niter=niter+1 !increase number of elapsed iterations by 1
    ENDDO
    
    IF(niter.EQ.itermax) THEN !if we stopped due to reaching maximum number of iterations, then stop
        PRINT*,' no convergence, increase nx or itermax in the f_xcce subroutine'
        STOP
    ENDIF
END FUNCTION

FUNCTION f_omega(surfg) result(omega) !relation between angular velocity and surface density for nearby (Kennicutt 1998) galaxies, as derived by Krumholz & McKee (2005)
    REAL, INTENT(in) :: surfg
    REAL :: msun,pc,Myr,surfg2,omega
    msun=1.989d30 !solar mass in kg
    pc=3.086d16 !parsec in meters
    Myr=1.e6*86400.*365.25 !million years in seconds
    surfg2=surfg/100./msun*pc**2. !surface density in units of 100 Msun/pc**2
    omega=0.058*surfg2**.49/Myr !angular velocity
END FUNCTION

FUNCTION f_cfe(surfg,qT,omega) result(cfearray) !FUNCTION TO CALCULATE THE CLUSTER FORMATION EFFICIENCY
    REAL, INTENT(in) :: surfg,qT,omega
    REAL :: pc,Msun,Myr
    REAL :: qvir,tsn,tview,surfGMC,ecore,beta0
    REAL :: mach,sigrho,mulnx,fbound,xsurv,fcce,fcce2
    REAL :: cfearray(1:4)
    INTEGER :: uparam,ioerror
    INTEGER :: sflaw,radfb
    CHARACTER*16 :: readfile
    
    !SET CONSTANTS
    pc=3.086e16 !parsec in meters
    Msun=1.989e30 !solar mass in kg
    Myr=1.e6*86400.*365.25 !million years in seconds
    
    !READ PARAMETERS
    uparam=18
    readfile='parameters.in'
    OPEN(UNIT=uparam,FILE=readfile,IOSTAT=ioerror)
    IF(ioerror.NE.0) THEN
        PRINT*,' stop -- cannot read parameters from file: ',readfile
        STOP
    ENDIF
    READ(uparam,*) sflaw !star formation law - NOTE: set to 0 for Elmegreen(2002) and to 1 for Krumholz & McKee (2005)
    READ(uparam,*) qvir !giant molecular cloud virial parameter
    READ(uparam,*) tsn !time of the first supernova
    READ(uparam,*) tview !time at which CFE is determined
    READ(uparam,*) surfGMC !giant molecular cloud surface density
    READ(uparam,*) ecore !maximum (protostellar core) star formation efficiency
    READ(uparam,*) beta0 !if turbulent-to-magnetic pressure ratio is not specified, set to turbulent-only
    READ(uparam,*) radfb !SN/radiative feedback mode - NOTE: set to 0 for supernovae only, to 1 for radiative only, and to 2 for both
    CLOSE(UNIT=uparam)
    tsn=tsn*Myr
    tview=tview*Myr
    surfGMC=surfGMC*Msun/pc**2.
    IF(surfg.GT.surfGMC) surfGMC=surfg

    !CALCULATE DERIVED PARAMETERS
    mach=f_mach(surfg,qT,omega) !Mach number
    sigrho=f_sigrho(mach,beta0) !dispersion of overdensity PDF
    mulnx=-.5*sigrho**2. !logarithmic mean of overdensity PDF

    !CALCULATE F_BOUND
    xsurv0=0.
    fbound=f_integrate(xsurv0,mulnx,sigrho,surfg,qT,omega,ecore,beta0,qvir,tsn,tview,surfGMC,0,sflaw,radfb) !naturally bound part of star formation

    !CALCULATE F_CCE
    xsurv=f_xcce(surfg,qT,omega,surfGMC,qvir,tview) !critical overdensity to remain bound despite the cruel cradle effect
    fcce=f_integrate(xsurv,mulnx,sigrho,surfg,qT,omega,ecore,beta0,qvir,tsn,tview,surfGMC,1,sflaw,radfb) !part of bound SF surviving the cruel cradle effect
    fcce2=f_integrate(xsurv,mulnx,sigrho,surfg,qT,omega,ecore,beta0,qvir,tsn,tview,surfGMC,2,sflaw,radfb) !part of all SF surviving the cruel cradle effect

    !CALCULATE CFE
    cfearray=[fbound*fcce,fbound,fcce,fcce2] !array containing the cluster formation efficiency, fbound, fcce, and fcce2 (i.e. fcce with respect to all SF)
END FUNCTION

END MODULE


