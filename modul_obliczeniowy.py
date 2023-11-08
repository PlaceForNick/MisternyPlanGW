import numpy as np
from math import *

a = 6378137.000
e2 = 0.00669438002290

class funkcje():
    
    def dms(self,x):
        '''
        
        Zmienia liczby w [rad] do formatu [° ' "] oraz wypluwa je na konsole.
    
        Parameters
        ----------
        x : liczba w [rad]
    
        Returns
        -------
        d : stopnie [°]
        m : minuty [']
        s : sekundy ["]
    
        '''
        znak = ' '
        if x < 0:
            znak = '-'
            x = abs(x)
        x = x * 180/pi
        d = int(x)
        m = int((x - d) *  60)
        s = (x - d - (m/60)) * 3600
        # print(znak, "%3d°%2d'%8.5f''" % (d, m, s))
        # return (d,m,s)
        return(f"{znak}{d:3d}\u00B0{m:2d}'{s:8.5f}''")
    
    def fromdms(self,X):
        '''
        
        Zmiana ze stopni w ukladzie [° ' "] na [rad] oraz [°]
    
        Parameters
        ----------
        X : liczba w formacie [° ' "]
            (zamiast ° ' " wstawiac spacje! Inaczej moze nie zadzialac)
    
        Returns
        -------
        Z : liczba w formacie [rad]
        Y : liczba w formacie [°]
    
        '''
        znak = 1
        if X[0] == '-':
             znak = -1
        Y = X.split(' ' or '"' or "'" or '°' or '-')
        d = int(Y[0])
        m = int(Y[1])
        s = float(Y[2])
        s = s/3600
        m = m/60
        Y = znak*(d+m+s)
        Z = Y * pi/180
        Y = float(f'{Y:7.5f}')
        return(Z,Y)
    
    def Np(self,f, a, e2):
        '''
        Liczy promien krzywizny w I wertykale (przekroj o MIN krzywiznie)
    
        Parameters
        ----------
        f : kat fi [rad]
        a : stala elipsoidy
        e2 : stala elipsoidy
    
        Returns
        -------
        N : promien krzywizny w I wertykale [m]
    
        '''
        N = a / np.sqrt(1 - e2 * np.sin(f)**2)
        return(N) 
    
    def Mp(self,f, a, e2):
        '''
        Liczy promien krzywizny w plaszczyznie poludnika (przekroj o MAX krzywiznie)
    
        Parameters
        ----------
        f : kat fi [rad]
        a : stala elipsoidy
        e2 : stala elipsoidy
    
        Returns
        -------
        M : promien krzywizny w I wertykale [m]
    
        '''
        M = ((a * (1 - e2)) / (np.sqrt((1 - e2 * (np.sin(f))**2)**3)))
        return(M)
    
    def xyz2flh(self,X, Y, Z, a, e2):
        '''
        [Algorytm HIRVONENA]
        
        Przelicza wspl XYZ na fi,lambda,h
    
        Parameters
        ----------
        X : wspl X [m]
        Y : wspl Y [m]
        Z : wspl Z [m]
        a : stala elipsoidy
        e2 : stala elipsoidy
    
        Returns
        -------
        f : wspl fi (szerokosc geo.) [rad]
        l : wspl lambda (dlugosc geo.) [rad]
        h : wspl h (wysokosc elipsoidalna) [m]
    
        '''
        p = np.sqrt(X**2 + Y**2) #promien rownoleznika
        f = np.arctan(Z/(p * (1 - e2))) #pierwsze przyblizenie f dla h=0
        while True:
            N = self.Np(f, a, e2)
            h = (p / np.cos(f)) - N
            fp = f
            f = np.arctan(Z / (p * (1 - e2 * N / (N + h))))
            if abs(fp - f) < (0.000001/206265):
                break
        
        N = self.Np(f,a,e2)
        h = (p / np.cos(f)) - N
        l = np.arctan2(Y, X)
        return(f, l, h)
    
    def flh2xyz(self,f,l,h,a,e2):
        '''
        Przelicza wspl fi,lambda,h na XYZ
    
        Parameters
        ----------
        f : wspl fi (szerokosc geo.) [rad]
        l : wspl lambda (dlugosc geo.) [rad]
        h : wspl h (wysokosc elipsoidalna) [m]
        a : stala elipsoidy
        e2 : stala elipsoidy
    
        Returns
        -------
        x : wspl X [m]
        y : wspl Y [m]
        z : wspl Z [m]
    
        '''
        N = self.Np(f,a,e2)
        x = (N+h)*np.cos(f)*np.cos(l)
        y = (N+h)*np.cos(f)*np.sin(l)
        z = ((N*(1-e2)+h))*np.sin(f)
        return(x,y,z) 
       
    def saz2neu(self,s, alfa, z):
        dneu = np.array([s * np.sin(z) * np.cos(alfa),
                         s * np.sin(z) * np.sin(alfa),
                         s * np.cos(z)])
        return(dneu)
    
    def neu2saz(self,dx):
        s = np.sqrt(dx @ dx)
        alfa = np.arctan2(dx[1],dx[2])
        z = np.arccos(dx[2]/s)
        return(s,alfa,z)
    
    def Rneu(self,f, l):
        R = np.array([[-np.sin(f) * np.cos(l), -np.sin(l), np.cos(f) * np.cos(l)],
                      [-np.sin(f) * np.sin(l), np.cos(l), np.cos(f) * np.sin(l)],
                      [np.cos(f), 0. ,np.sin(f)]])
        return(R)
        
    
    def neu2xyz(self,dneu, f, l):
        R = Rneu(f, l)
        dXYZ = R @ dneu
        return(dXYZ)
    
    def xyz2neu(self,dX,f,l):
        R = Rneu(f, l)
        return(R.T @ dX)
        
    def kivioj(self,f,l,A,s,a,e2): #wspolrzednje PUNKTU A
        '''
        Sluzy do liczenia zadania wprost przy uzyciu algorytmu kivioja.
        Przelicza wspl pkt poczatkowego na wspl pkt koncowego

        Parameters
        ----------
        f : wspl fi pkt poczatkowego (szerokosc geo.) [rad]
        l : wspl lambda pkt koncowego  (dlugosc geo.) [rad]
        A : azymut wprost (poczatek-koniec) [rad]
        s : dl lini geodezyjnej (odleglosc elipsoidalna) [m]
        a : stala elipsoidy
        e2 : stala elipsoidy

        Returns
        -------
        f : wspl fi pkt koncowego(szerokosc geo.) [rad]
        l : wspl lambda pkt koncowego (dlugosc geo.) [rad]
        A : azymut odwrotny (koniec-poczatek) [rad]

        '''
        n = int(s/1000)
        ds = s/n
        c1 = self.Np(f,a,e2) * np.cos(f) * np.sin(A)
        for i in range(n):
            M = self.Mp(f,a,e2)
            N = self.Np(f,a,e2)
            df = ds * np.cos(A)/M
            dA = ds * np.sin(A) * np.tan(f) / N
            fm = f + df/2
            Am = A + dA/2
            Mm = self.Mp(fm,a,e2)
            Nm = self.Np(fm,a,e2)
            dfm = (ds * np.cos(Am))/Mm
            dAm = (ds * np.sin(Am) * np.tan(fm))/Nm
            dlm = (ds * np.sin(Am))/(Nm * np.cos(fm))
            f = f + dfm
            l = l + dlm
            A = A + dAm
            c2 = self.Np(f,a,e2) * np.cos(f) * np.sin(A)
        A = A + pi
        if A > 2*pi:
            A = A - 2*pi
        return(f,l,A,c1,c2)  #wspolrzednje PUNKTU B
    
    def vincenty(self,fa,la,fb,lb,a,e2):
        '''
        Sluzy do liczenia zadania odwrotnego przy uzyciu algorytmu vincentego.
        Przelicza wspl pkt poczatkowego i koncowego na azymut i linie geodezyjna.


        Parameters
        ----------
        fa : wspl fi pkt poczatkowego (szerokosc geo.) [rad]
        la : wspl lambda pkt koncowego  (dlugosc geo.) [rad]
        fb : wspl fi pkt koncowego (szerokosc geo.) [rad]
        lb : wspl lambda pkt koncowego (dlugosc geo.) [rad]
        a : stala elipsoidy
        e2 : stala elipsoidy

        Returns
        -------
        s : dl lini geodezyjnej (odleglosc elipsoidalna) [m]
        Aab : azymut wprost (poczatek-koniec) [rad]
        Aba : azymut odwrotny (koniec-poczatek) [rad]

        '''
        b = a * np.sqrt(1 - e2)
        fl = 1 - (b / a) #splaszczenie
        Ua = np.arctan((1 - fl) * np.tan(fa)) #szer zredukowana pkt A
        Ub = np.arctan((1 - fl) * np.tan(fb)) #szer zredukowana pkt B
        dl = lb - la
        L = dl
        while True:
            sin_sigma = np.sqrt((np.cos(Ub) * np.sin(L))**2 + ((np.cos(Ua) * np.sin(Ub)) - (np.sin(Ua) * np.cos(Ub) * np.cos(L)))**2)
            cos_sigma = (np.sin(Ua) * np.sin(Ub)) + (np.cos(Ua) * np.cos(Ub) * np.cos(L))
            sigma = np.arctan(sin_sigma/cos_sigma)
            sin_alfa = (np.cos(Ua) * np.cos(Ub) * np.sin(L)) / sin_sigma
            cos2_alfa = 1 - sin_alfa**2
            cos2_sigma_m = cos_sigma - (2*(np.sin(Ua) * np.sin(Ub)) / cos2_alfa)
            C = (fl/16) * cos2_alfa * (4 + fl*(4 - 3 * cos2_alfa)) 
            Ls = L
            L = dl + (1-C) * fl * sin_alfa * (sigma + C * sin_sigma * (cos2_sigma_m + C * cos_sigma * (-1 + 2*(cos2_sigma_m**2))))
            if abs(Ls - L) < (0.000001/206265):
                break
        u2 = ((a**2 - b**2) / (b**2)) * cos2_alfa
        A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175*u2)))
        B = (u2/1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
        d_sigma = B * sin_sigma * (cos2_sigma_m + ((1/4) * B * ((cos_sigma * (-1+2*cos2_sigma_m**2) - ((1/6) * B * cos2_sigma_m * (-3 + 4*sin_sigma**2) * (-3 + 4 * cos2_sigma_m**2))))))   
        sab = b * A * (sigma - d_sigma)
        Aab = np.arctan2((np.cos(Ub) * np.sin(L)) , (( np.cos(Ua) * np.sin(Ub)) - (np.sin(Ua) * np.cos(Ub) * np.cos(L))))
        Aba = np.arctan2((np.cos(Ua) * np.sin(L)) , ((-np.sin(Ua) * np.cos(Ub)) + (np.cos(Ua) * np.sin(Ub) * np.cos(L)))) + pi
        if Aab > 2*pi:
            Aba = Aba - 2*pi
        if Aba < 0:
            Aba = Aba + 2*pi
        if Aab < 0:
            Aab = Aab + 2*pi
        return(sab,Aab,Aba)
    
    def sigma(self,f, a, e2):
        A0 = 1 - e2/4 - 3 * e2**2/64 - 5 * e2**3/256
        A2 = (3/8) * (e2 + e2**2/4 + 15*e2**3/128)
        A4 = (15/256) * (e2**2 + (3 * e2**3)/4)
        A6 = 35 * e2**3/3072
        sigma = a * (A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))
        return(sigma)
    
    def f1(self,xgk, a, e2):
        A0 = 1 - e2/4 - 3 * e2**2/64 - 5 * e2**3/256
        f = xgk / (a * A0)
        while True:
            fs = f
            s = self.sigma(f, a, e2)
            f = fs + ((xgk - s)/(a * A0))
            if abs(fs-f) < (0.000001/206265):
                break
        return(f)
    
    def fl2PL1992(self,f,l,a,e2,l0=radians(19), m0 = 0.9993):
        b2 = a**2*(1 - e2)
        ep2 = (a**2 - b2)/b2
        dl = l - l0
        t = tan(f)
        n2 = ep2 * cos(f)**2
        N = self.Np(f,a,e2)
        sigm = self.sigma(f,a,e2)
        xgk = sigm + (dl**2/2) * N * sin(f)*cos(f)*(1 + (dl**2/12)*cos(f)**2*(5-t**2+9*n2+4*n2**2)+ ((dl**4)/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
        ygk = dl*N*cos(f)*(1+(dl**2/6)*cos(f)**2*(1 - t**2 + n2) + (dl**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
        x92 = xgk * m0 - 5300000
        y92 = ygk * m0 + 500000
        return(x92,y92,xgk,ygk)
        
    def fl2GK(self,fa,la,a,e2,l0=radians(19)):
        b2=a**2*(1-e2)
        e_2=(a**2-b2)/b2
        dl=la-l0
        t=tan(fa)
        n2=e_2*cos(fa)**2
        N=Np(fa,a,e2)
        sigm=sigma(fa,a, e2)
        x=sigm+(dl**2/2)*N*sin(fa)*cos(fa)*(1+(dl**2/12)*(cos(fa))**2*(5-t**2+9*n2+4*n2**2)+(dl**4/360)*(cos(fa))**4*(61-58*t**2+t**4+270*n2-330*n2*t**2))
        y=dl*N*cos(fa)*(1+(dl**2/6)*(cos(fa))**2*(1-t**2+n2)+(dl**4/120)*(cos(fa))**4*(5-18*t**2+t**4+14*n2-58*n2*t**2))
        return(x,y)
    
    def PL19922fl(self,x,y,a,e2,l0=radians(19),m0=0.9993):
        xgk = (x + 5300000)/m0
        ygk = (y - 500000)/m0
        fi1 = self.f1(xgk,a,e2)
        N1 = self.Np(fi1,a,e2)
        M1 = self.Mp(fi1,a,e2)
        t = tan(fi1)
        b2 = a**2*(1-e2)
        ep2 = (a**2 - b2)/b2
        n2 = ep2 * cos(fi1)**2
        fi = fi1 - ((ygk**2*t)/(2*M1*N1)) * (1 - (ygk**2/(12*N1**2)) * (5 + 3*t**2 + n2 - 9*n2*t**2 - 4*n2**20) + (ygk**4/(360*N1**4))*(61 + 90*t**2 + 45*t**4))
        lam = l0 + (ygk/(N1*cos(fi1)))*(1 - (ygk**2/(6*N1**2))*(1 + 2*t**2 + n2) + (ygk**4/(120*N1**4)) * (5 + 28*t**2 + 24*t**4 +6*n2 + 8*n2*t**2))
        return(fi,lam,xgk,ygk)
    
    def GK2fl(self,x,y,a,e2,l0=radians(19)):
        b2=a**2*(1-e2)
        e_2=(a**2-b2)/b2
        A0=1-e2/4-3*e2**2/64-5*e2**3/256
        fl=x/(a*A0)
        while True:
            fs=fl
            sigm=sigma(fl,a,e2)
            fl=fl+(x-sigm)/(a*A0)
            if abs( fl-fs)<(0.000001/206265):
                break
        N=Np(fl,a,e2)
        M=Mp(fl,a,e2)
        t=tan(fl)
        n2=e_2*(cos(fl)**2)
        f= fl - (((y**2) * t)/(2 * M * N)) * (1 - ((y**2)/(12*N**2)) * (5 + 3 * t**2 + n2 - 9 * n2 * t**2 - 4 * n2**2) + ((y**4)/(360 * N**4)) * (61 + 90 * t**2 + 45 * t**4))
        l=l0+(y/(N*cos(fl)))*(1-((y**2)/(6*N**2))*(1+2*t**2+n2)+((y**4)/(120*N**4))*(5+28*t**2+24*t**4+6*n2+8*n2*t**2))
        return(f,l)
       
    def fl2PL2000(self,f,l,a,e2,ns,m0= 0.999923):
        if ns == 5:
            l0 = radians(15)
        elif ns == 6:
            l0 = radians(18)
        elif ns == 7:
            l0 = radians(21)
        elif ns == 8:
            l0 = radians(24)
        b2 = a**2*(1 - e2)
        ep2 = (a**2 - b2)/b2
        dl = l - l0
        t = tan(f)
        n2 = ep2 * cos(f)**2
        N = self.Np(f,a,e2)
        sigm = self.sigma(f,a,e2)
        xgk = sigm + (dl**2/2) * N * sin(f)*cos(f)*(1 + (dl**2/12)*cos(f)**2*(5-t**2+9*n2+4*n2**2)+ ((dl**4)/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
        ygk = dl*N*cos(f)*(1+(dl**2/6)*cos(f)**2*(1 - t**2 + n2) + (dl**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
        x2000 = xgk * m0
        y2000 = ygk * m0 + ns * 1000000 + 500000
        return(x2000,y2000,xgk,ygk)
    
    def strefa(self,l):
        '''
        Automatycznie wybiera strefe odwzorowawcza dla ukladu wspolrzednych PL2000

        Parameters
        ----------
        l : wspl lambda (dlugosc geo.) [rad]

        Returns
        -------
        ns : nr strefy [-]

        '''
        if l < radians(16.5) and l > radians(13.5): #ns = 5
            l0 = radians(15)
            ns = 5
        elif l < radians(19.5) and l > radians(16.5): #ns = 6
            l0 = radians(18)
            ns = 6
        elif l < radians(22.5) and l > radians(19.5): #ns = 7
            l0 = radians(21)
            ns = 7
        elif l < radians(25.5) and l > radians(22.5): #ns = 8
            l0 = radians(24)
            ns = 8
        else:
            raise NotImplementedError('Wyznaczona strefa jest nieprawidlowa dla odwzorowania PL2000.'
                                      'Sprawdz poprawnosc wspolrzednej l. Podana wartosc to:'
                                      f'l = {self.dms(l)}')
        return(ns)
    
    def strefa2(self,y2000):
        '''
        Automatycznie wybiera strefe odwzorowawcza dla ukladu wspolrzednych PL2000

        Parameters
        ----------
        y2000 : wspl y w ukladdzie PL2000 [m]

        Returns
        -------
        ns : nr strefy [-]

        '''
        ns = str(y2000)[0]
        if ns != '5' and '6' and '7' and '8':
            raise NotImplementedError('Wyznaczona strefa jest nieprawidlowa dla odwzorowania PL2000.'
                                      'Sprawdz poprawnosc wspolrzednej y2000. Podana wartosc to:'
                                      f'y2000 = {y2000:0.3f}')
        return(int(ns))
            
    def PL20002fl(self,x20,y20,a,e2,ns,m0 = 0.999923):
        if ns == 5:
            l0 = radians(15)
        elif ns == 6:
            l0 = radians(18)
        elif ns == 7:
            l0 = radians(21)
        elif ns == 8:
            l0 = radians(24)
        xgk = x20/m0
        ygk = (y20 - 500000 - ns*1000000)/m0
        fi1 = self.f1(xgk,a,e2)
        N1 = self.Np(fi1,a,e2)
        M1 = self.Mp(fi1,a,e2)
        t = tan(fi1)
        b2 = a**2*(1-e2)
        ep2 = (a**2 - b2)/b2
        n2 = ep2 * cos(fi1)**2
        fi = fi1 - ((ygk**2*t)/(2*M1*N1)) * (1 - (ygk**2/(12*N1**2)) * (5 + 3*t**2 + n2 - 9*n2*t**2 - 4*n2**20) + (ygk**4/(360*N1**4))*(61 + 90*t**2 + 45*t**4))
        lam = l0 + (ygk/(N1*cos(fi1)))*(1 - (ygk**2/(6*N1**2))*(1 + 2*t**2 + n2) + (ygk**4/(120*N1**4)) * (5 + 28*t**2 + 24*t**4 +6*n2 + 8*n2*t**2))
        return(fi,lam,xgk,ygk)
    
    def mgk(self,xgk,ygk,a,e2):
        f = f1(xgk,a,e2)
        R = np.sqrt(Np(f, a, e2) * Mp(f, a, e2))
        m = 1 + ygk**2 / (2 * R**2) + ygk**4 / (24 * R**4)
        return(m)
    
    # def red_gk(xa,ya,xb,yb,a,e2,l0 = radians(19)):
    #     xm = (xa + xb)/2
    #     ym = (ya + yb)/2
    #     fm,lm = GK2fl(xm,ym,a,e2,l0)
    #     Rm2 = Mp(fm,a,e2) * Np(fm,a,e2)
    #     print(np.sqrt(Rm2))
    #     sgk = np.sqrt((xa-xb)**2 + (ya-yb)**2)
    #     r = sgk * (ya**2 + ya*yb + yb**2)/(6 * Rm2)
    #     selip = sgk - r
    #     return(r,selip,sgk)
    
    def red_gk(self,xa,ya,xb,yb,ns,a,e2):
        xm = (xa + xb)/2
        ym = (ya + yb)/2
        if ns == 'GK':
            l0 = 0
            fm,lm = GK2fl(xm, ym, a, e2, l0)
        elif ns =='1992':
            fm, lm, xm, ym = PL19922fl(xm, ym, a, e2)
            fa, la, xa, ya = PL19922fl(xa, ya, a, e2)
            fb, lb, xb, yb = PL19922fl(xb, yb, a, e2)
        elif (ns == 5) or (ns == 6) or (ns == 7) or (ns == 8):
            fm, lm, xm, ym = PL20002fl(xm, ym, a, e2, ns)
            fa, la, xa, ya = PL20002fl(xa, ya, a, e2, ns)
            fb, lb, xb, yb = PL20002fl(xb, yb, a, e2, ns)
        else:
            fm,lm = GK2fl(xm, ym, a, e2, ns)
        Rm2 = Mp(fm,a,e2)*Np(fm,a,e2)
        sgk = np.sqrt((xa - xb)**2 + (ya - yb)**2)
        r = sgk *(ya**2 + ya*yb + yb**2)/(6 * Rm2)
        selip = sgk - r 
        return(r,selip,sgk)
    
    def redu_d(self,xa, xb, ya, yb, ha, hb, spom, l0, a, e2):
        xm = (xa + xb) / 2
        ym = (ya + yb) /2
        fm, lm = GK2fl(xm,ym,a,e2,l0)
        Rm = sqrt((Np(fm,a,e2))*(Mp(fm,a,e2)))
        dh = hb-ha
        s0 = sqrt((spom**2 - dh**2)/((1+(ha/Rm))*(1+(hb/Rm))))
        selip = 2 * Rm * asin(s0/(2*Rm))
        return selip
    
    def red_odl_skosnej(self,x,y,h,s,z,beta,i,l,a,e2):
        Xp = x + s*np.sin(z) *np.cos(beta)
        Yp = y + s*np.sin(z) *np.sin(beta)
        hp = h + i + s*np.cos(z) - l
        Xm = (x + Xp)/2
        Ym = (y + Yp)/2
        fm, lm = GK2fl(Xm,Ym,a,e2)
        Nm = Np(fm, a, e2)
        Mm = Mp(fm, a, e2)
        Rm = np.sqrt(Nm * Mm)
        Ha = h  + i
        Hb = hp + l
        s0 = np.sqrt((s**2 - (Hb - Ha)**2)/(1+(Ha/Rm)*(1 + (Hb/Rm))))
        return(s0, Xp, Yp)
    
    def zbiez_pld(self,xgk,ygk,a,e2):
        fi1 = f1(xgk,a,e2)
        N1 = Np(fi1,a,e2)
        t = tan(fi1)
        b2 = a**2*(1-e2)
        ep2 = (a**2 - b2)/b2
        n2 = ep2 * (cos(fi1))**2
        gamma = (ygk / N1) * t * (1 - (ygk**2 / (3 * N1**2)) * (1 + t**2 - n2 - 2 * n2**2) + (ygk**4 / (15 * N1**4)) * (2 + 5 * t**2 + 3 * t**4))
        return(gamma)
    
    # def zbiez_pld(xgk,ygk,a,e2):
    #     fii = f1(xgk, a, e2)
    #     N1 = Np(fii, a, e2)
    #     b2 = a**2 * (1-e2) 
    #     e2p = ( a**2 - b2 ) / b2 
    #     t1 = np.tan(fii)
    #     ni = np.sqrt(e2p * (np.cos(fii))**2)
    #     #g = ygk/N1*t1*(1 - (ygk)**2/3*(N1**2)*(1 + t1**2 - ni**2 -2*ni**2) + (ygk**4/15*N1**4)*(2 + 5*t1**2 + 3*t1**4))
    #     g = ygk/Np(fii,a,e2) * t1 * (1 - ygk**2/(3*Np(fii,a,e2)**2) * (1 + t1**2 - ni - 2 * ni**2) + ygk**4/(15*Np(fii,a,e2)**4 *(2 + 5 * t1**2 + 3 * t1**4)))
    #     return(g) #gamma
    
    def zbiez_kier(self,xa,ya,xb,yb,a,e2,l0):
        xm = (xa + xb)/2
        ym = (ya + yb)/2
        fm,lm = GK2fl(xm,ym,a,e2,l0)
        M = Mp(fm,a,e2) 
        N = Np(fm,a,e2)
        Rm2 = M*N
        delta = ((xb - xa) * (2 * ya + yb)) / (6 * Rm2)
        # dab = ((xb - xa)*(2*ya + yb))/(6*Rm2)
        # print(dab)
        # dba = ((xa - xb)*(2*yb + ya))/(6 *Rm2)
        return(delta)
    
    # def zbiez_kier(xa,ya,xb,yb,l0,a,e2):
    #     x = (xa + xb)/2
    #     y = (ya + yb)/2
    #     f,l = GK2fl(x, y, l0, a, e2)
    #     M = Mp(f,a,e2) 
    #     N = Np(f,a,e2)
    #     R = np.sqrt(M*N)
    #     dab = (xb - xa)*(2*ya + yb)/(6*R**2)
    #     dba = (xa - xb)*(2*yb + ya)/(6 *R**2)
    #     return(dab,dba) #deltaAB deltaBA
    
    def azymut(self,xa,ya,xb,yb):
        alfa_ab = np.arctan(((yb - ya) / (xb - xa)))
        alfa_ba = alfa_ab + np.pi
        if alfa_ab < 0:
            alfa_ab += 2*pi
        if alfa_ba < 0:
            alfa_ba += 2*pi
        return(alfa_ab,alfa_ba)
    
    def redukcja_AZ(self,xa,ya,xb,yb,dAb,dBa,a,e2):
        # dAb = zbiez_kier(xa,ya,xb,yb,l0,a,e2)
        # dBa = zbiez_kier(xb,yb,xa,ya,l0,a,e2)
        alfa_ab, alfa_ba = azymut(xa,ya,xb,yb)
        AAb = alfa_ab + zbiez_pld(xa, ya, a, e2) + dAb
        ABa = alfa_ba + zbiez_pld(xb, yb, a, e2) + dBa
        if alfa_ab < 0:
            alfa_ab += 2*pi
        if alfa_ba < 0:
            alfa_ba += 2*pi
        if AAb < 0:
            AAb += 2*pi
        if ABa < 0:
            ABa += 2*pi
        return(alfa_ab,alfa_ba,AAb,ABa)
    
    def bursa(self,x,p):
        '''
        transformacja Bursa-Wolfa / konforenna / przez podobienstwo
        - znieksztalcony orot, skala
        - zachowany ksztalt figury (katy)
        
        Parameters
        ----------
        p : [kappa, alfa, beta, gamma, dx, dy, dz]
        x : [x, y, z]
    
        Returns
        -------
        Xw : [xw, yw, zw] - wspolrzedne wtorne
    
        '''
        A = np.array([[ p[0],  p[3], -p[2]],
                      [-p[3],  p[0],  p[1]],
                      [ p[2], -p[1],  p[0]]])
        T = np.array([p[4], p[5], p[6]])
        Xw = x + A@x + T
        return(Xw)
    
    def kwazi(self,x,p):
        '''
        transformacja kwaziafiniczna
        -znieksztalcony obrot, skala, ksztalt figury (katy)
    
        Parameters
        ----------
        p : [kappa x, kapppa y, kappa z, alfa, beta, gamma, dx, dy, dz]
        x : [x, y, z]
    
        Returns
        -------
        Xw : [xw, yw, zw] - wspolrzedne wtorne
    
        '''
        A = np.array([[ p[0],  p[5], -p[4]],
                      [-p[5],  p[1],  p[3]],
                      [ p[4], -p[3],  p[2]]])
        T = np.array([p[6], p[7], p[8]])
        Xw = x + A@x + T
        return(Xw)
    
    def izometr(self,x,p):
        '''
        transformacja izometryczna / sztywnego ciala
        - zachowany ksztalt figury (katy), skala
        - znieksztalcony obrot
        
        Parameters
        ----------
        p : [alfa, beta, gamma, dx, dy, dz]
        x : [x, y, z]
    
        Returns
        -------   
        Xw : [xw, yw, zw] - wspolrzedne wtorne
    
        '''
        A = np.array([[ 0   ,  p[2], -p[1]],
                      [-p[2],  0   ,  p[0]],
                      [ p[1], -p[0],  0   ]])
        T = np.array([p[3], p[4], p[5]])
        Xw = x + A@x + T
        return(Xw)
    
