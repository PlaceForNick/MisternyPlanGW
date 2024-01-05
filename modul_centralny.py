import numpy as np
from math import *
from modul_obliczeniowy import funkcje

class Format:
    normal = '\033[0m'
    podkresl = '\033[4m'
    kursywa = '\x1B[3m'
    
class skrypt(funkcje):
    
    ID = 0
    
    def __init__(self, model='grs80', zapis=False, posrednie=False, nazwa='', X='', Y='', Z='', f='', l='', h='', X2='', Y2='', Z2='', s='', A='', z ='', x2000='', y2000='', x1992='', y1992='', xgk='', ygk='', f2='', l2='', h2='', s_elip='', x2000_2='', y2000_2='', x1992_2='', y1992_2='', xgk2='', ygk2='', s_0='', s_1992='', s_2000='', p=''):

        skrypt.ID += 1
        
        self.__elipsoida(model) #wybor elipsoidy
        # self.__zapiszplik(zapis, nazwa) #wybor zapisu do pliku txt     
        self.posrednie = posrednie #definiuje czy ma wypluwac obliczenia posrednie
        
        self.m0_1992 = 0.9993
        self.m0_2000 = 0.999923
        
        #zamiana podanych danych na liste
        dane = [X, Y, Z, f, l, h, X2, Y2, Z2, s, A, z, x2000, y2000, x1992, y1992, xgk, ygk, f2, l2, h2, s_elip, x2000_2, y2000_2, x1992_2, y1992_2, xgk2, ygk2, s_0, s_1992, s_2000, p]
        dane_ost = []
        for wartosc in dane:
            if wartosc == p and p!= '':
                wartosc_lista = []
                wartosc_lista.append(wartosc)
            elif type(wartosc) == list:
                wartosc_lista = wartosc
            else:
                wartosc_lista = []
                wartosc_lista.append(wartosc)
            dane_ost.append(wartosc_lista)
        self.X = dane_ost[0]
        self.Y = dane_ost[1]
        self.Z = dane_ost[2]
        self.f = dane_ost[3]
        self.l = dane_ost[4]
        self.h = dane_ost[5]
        self.X2 = dane_ost[6]
        self.Y2 = dane_ost[7]
        self.Z2 = dane_ost[8]
        self.s = dane_ost[9]
        self.A = dane_ost[10]
        self.z = dane_ost[11]
        self.x2000 = dane_ost[12]
        self.y2000 = dane_ost[13]
        self.x1992 = dane_ost[14]
        self.y1992 = dane_ost[15]
        self.xgk = dane_ost[16]
        self.ygk = dane_ost[17]
        self.f2 = dane_ost[18]
        self.l2 = dane_ost[19]
        self.h2 = dane_ost[20]
        self.s_elip = dane_ost[21]
        self.x2000_2 = dane_ost[22]
        self.y2000_2 = dane_ost[23]
        self.x1992_2 = dane_ost[24]
        self.y1992_2 = dane_ost[25]
        self.xgk2 = dane_ost[26]
        self.ygk2 = dane_ost[27]
        self.s_0 = dane_ost[28]
        self.s_1992 = dane_ost[29]
        self.s_2000 = dane_ost[30]
        self.p = dane_ost[31]
        
        dane_kat = [self.f, self.l, self.A, self.z, self.f2, self.l2, self.p]
        dane_kat_ost = []
        for wartosc_lista in dane_kat:
            wartosc_ost = []
            i = 0  
            while True:
                try:
                    wartosc = wartosc_lista[i]
                    if type(wartosc) == list and len(wartosc) == 1:
                        wartosc = wartosc[0]
                    if type(wartosc) == list:
                        if len(wartosc) == 9:
                            dane_kat.append([wartosc[3]])
                            dane_kat.append([wartosc[4]])
                            dane_kat.append([wartosc[5]])
                        elif len(wartosc) == 7:
                            dane_kat.append([wartosc[1]])
                            dane_kat.append([wartosc[2]])
                            dane_kat.append([wartosc[3]])
                        elif len(wartosc) == 6:
                            dane_kat.append([wartosc[0]])
                            dane_kat.append([wartosc[1]])
                            dane_kat.append([wartosc[2]])
                    elif wartosc =='':
                        wartosc = wartosc
                    elif type(wartosc) == str:
                        wartosc = self.fromdms(wartosc)[0]
                    else:
                        wartosc = radians(wartosc)
                    wartosc_ost.append(wartosc)
                    i += 1
                except IndexError:
                    break
            dane_kat_ost.append(wartosc_ost)
        self.f = dane_kat_ost[0]
        self.l = dane_kat_ost[1]
        self.A = dane_kat_ost[2]
        self.z = dane_kat_ost[3]
        self.f2 = dane_kat_ost[4]
        self.l2 = dane_kat_ost[5]
        
        p_nowe = []
        for parametry in self.p:
            if len(parametry) == 9:
                parametry[3:6] = dane_kat_ost[7][0], dane_kat_ost[8][0], dane_kat_ost[9][0]
            elif len(parametry) == 7:
                parametry[1:4] = dane_kat_ost[7][0], dane_kat_ost[8][0], dane_kat_ost[9][0]
            elif len(parametry) == 6:
                parametry[0:3] = dane_kat_ost[7][0], dane_kat_ost[8][0], dane_kat_ost[9][0]
            p_nowe.append(parametry)
        self.p = p_nowe
        
    def __elipsoida(self, model):
        #wybor elipsoidy
        if    model  == 'kra':
            self.a= 6378245
            self.b= 6356863.01877
        elif  model == "wgs84":
            self.a = 6378137.0 
            self.b = 6356752.31424518 
        elif  model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            # raise NieprawidlowaWartosc(f"{model} ten model elipsoidy nie jest obslugiwany")
            pass
        self.splasz = (self.a - self.b) / self.a
        self.e2 = (2 * self.splasz - self.splasz ** 2)
        
    def flh(self):
        '''
        Przelicza do ukladu wspolrzednych krzywoliniowych flh na podstawie dostepnych danych:
            -XYZ
            -PL2000
            -PL1992

        Returns
        -------
        f [deg]
        l [deg]
        h [m]

        '''   
        f_st = []; l_st = []
        f_ost = []; l_ost = []
        h_ost = []
        i = 0   
        
        if self.X != [''] and self.Y != [''] and self.Z != ['']:
            
            while i < len(self.X):
                
                f,l,h = self.xyz2flh(self.X[i], self.Y[i], self.Z[i], self.a, self.e2)
                f_st.append(self.dms(f)); l_st.append(self.dms(l))
                f_ost.append(np.rad2deg(f)); l_ost.append(np.rad2deg(l))
                h_ost.append(h)
                i += 1
        
        elif self.x2000 != [''] and self.y2000 != ['']:
            
            while i < len(self.x2000):

                ns = int(str(self.y2000[i])[0])
                f,l,xgk,ygk = self.PL20002fl(self.x2000[i], self.y2000[i], self.a, self.e2, ns)
                f_st.append(self.dms(f)); l_st.append(self.dms(l))
                f_ost.append(np.rad2deg(f)); l_ost.append(np.rad2deg(l))
                h_ost = self.h
                i += 1
        
        elif self.x1992 != [''] and self.y1992 != ['']:
            
            while i < len(self.x1992):

                f,l,xgk,ygk = self.PL19922fl(self.x1992[i], self.y1992[i], self.a, self.e2)
                f_st.append(self.dms(f)); l_st.append(self.dms(l))
                f_ost.append(np.rad2deg(f)); l_ost.append(np.rad2deg(l))
                h_ost = self.h
                i += 1
        
        if h_ost != ['']:
            print(Format.podkresl + f'\nZapytanie {skrypt.ID} [flh]:' + Format.normal + '\n\u03C6: ',('{} '*len(f_st)).format(*f_st), '\n\u03BB: ',('{} '*len(l_st)).format(*l_st), '\nh: ',('{:.3f} '*len(h_ost)).format(*h_ost), '[m]')
        else:
            print(Format.podkresl + f'\nZapytanie {skrypt.ID} [flh]:' + Format.normal + '\n\u03C6: ',('{} '*len(f_st)).format(*f_st), '\n\u03BB: ',('{} '*len(l_st)).format(*l_st))
 
        return(f_ost,l_ost,h_ost)
    
    def XYZ(self):
        '''
        Przelicza do ukladu wspolrzednych prostokatnych XYZ na podstawie dostepnych danych:
            -flh
            -PL2000
            -PL1992

        Returns
        -------
        X [m]
        Y [m]
        Z [m]

        '''
        X_ost = []; Y_ost = []; Z_ost = []
        i = 0
        
        if self.f != [''] and self.l != [''] and self.h != ['']:
            
            while i < len(self.f):
                
                X,Y,Z = self.flh2xyz(self.f[i], self.l[i], self.h[i], self.a, self.e2)
                X_ost.append(X); Y_ost.append(Y); Z_ost.append(Z)
                i += 1
         
        elif self.x2000 != [''] and self.y2000 != [''] and self.h != ['']:
            
            while i < len(self.x2000):
                
                ns = int(str(self.y2000[i])[0])
                f,l,xgk,ygk = self.PL20002fl(self.x2000[i], self.y2000[i], self.a, self.e2, ns)
                X,Y,Z = self.flh2xyz(f, l, self.h[i], self.a, self.e2)
                X_ost.append(X); Y_ost.append(Y); Z_ost.append(Z)
                i += 1
                
        elif self.x1992 != [''] and self.y1992 != [''] and self.h != ['']:
            
            while i < len(self.x1992):
                
                f,l,xgk,ygk = self.PL19922fl(self.x1992[i], self.y1992[i], self.a, self.e2)
                X,Y,Z = self.flh2xyz(f, l, self.h[i], self.a, self.e2)
                X_ost.append(X); Y_ost.append(Y); Z_ost.append(Z)
                i += 1
                
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [XYZ]:' + Format.normal + '\nX: ',('{:.3f} '*len(X_ost)).format(*X_ost), '[m]\nY: ',('{:.3f} '*len(Y_ost)).format(*Y_ost), '[m]\nZ: ',('{:.3f} '*len(Z_ost)).format(*Z_ost), '[m]' )
        return(X_ost,Y_ost,Z_ost)
    
    def PL2000(self):
        '''
        Przelicza do ukladu wspolrzednych PL2000 na podstawie dostepnych danych:
            -XYZ
            -fl
            -PL1992

        Returns
        -------
        x2000 [m]
        y2000 [m]

        '''
        x_ost = []; y_ost = []
        i = 0
        
        if self.f != [''] and self.l != ['']:
            print(self.f)
            print(self.l)
            while i < len(self.f):
                
                x2000,y2000,xgk,ygk = self.fl2PL2000(self.f[i], self.l[i], self.a, self.e2, self.strefa(self.l[i]))
                x_ost.append(x2000); y_ost.append(y2000)
                i += 1
                
        elif self.X != [''] and self.Y != [''] and self.Z != ['']:
            
            while i < len(self.X):
                
                f,l,h = self.xyz2flh(self.X[i], self.Y[i], self.Z[i], self.a, self.e2)
                x2000,y2000,xgk,ygk = self.fl2PL2000(f, l, self.a, self.e2, self.strefa(l))
                x_ost.append(x2000); y_ost.append(y2000)
                i += 1
                
        elif self.x1992 != [''] and self.y1992 != ['']:
            
            while i < len(self.x1992):
                
                f,l,xgk,ygk = self.PL19922fl(self.x1992[i], self.y1992[i], self.a, self.e2)
                x2000,y2000,xgk,ygk = self.fl2PL2000(f, l, self.a, self.e2, self.strefa(l))
                x_ost.append(x2000); y_ost.append(y2000)
                i += 1
                
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [PL2000]:' + Format.normal + '\nx2000: ',('{:.3f} '*len(x_ost)).format(*x_ost), '[m]\ny2000: ',('{:.3f} '*len(y_ost)).format(*y_ost), '[m]' )
        return(x_ost,y_ost)
    
    def PL1992(self):
        '''
        Przelicza do ukladu wspolrzednych PL2000 na podstawie dostepnych danych:
            -XYZ
            -fl
            -PL2000

        Returns
        -------
        x1992 [m]
        y1992 [m]

        '''
        x_ost = []; y_ost = []
        i = 0
        
        if self.f != [''] and self.l != ['']:
            
            while i < len(self.f):
                
                x1992,y1992,xgk,ygk = self.fl2PL1992(self.f[i], self.l[i], self.a, self.e2)
                x_ost.append(x1992); y_ost.append(y1992)
                i += 1
                
        elif self.X != [''] and self.Y != [''] and self.Z != ['']:
            
            while i < len(self.X):
                
                f,l,h = self.xyz2flh(self.X[i], self.Y[i], self.Z[i], self.a, self.e2)
                x1992,y1992,xgk,ygk = self.fl2PL1992(f, l, self.a, self.e2)
                x_ost.append(x1992); y_ost.append(y1992)
                i += 1
                
        elif self.x2000 != [''] and self.y2000 != ['']:
            
            while i < len(self.x2000):
                
                f,l,xgk,ygk = self.PL20002fl(self.x2000[i], self.y2000[i], self.a, self.e2, self.strefa2(self.y2000[i]))
                x1992,y1992,xgk,ygk = self.fl2PL1992(f, l, self.a, self.e2)
                x_ost.append(x1992); y_ost.append(y1992)
                i += 1    
                
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [PL1992]:' + Format.normal + '\nx1992: ',('{:.3f} '*len(x_ost)).format(*x_ost), '[m]\ny1992: ',('{:.3f} '*len(y_ost)).format(*y_ost), '[m]' )
        return(x_ost,y_ost)
        
    def wprost(self):
        '''
        Wykonuje zadanie wprost przy uzyciu ALGORYTMU KIVIOJA na podstawie dostepnych danych:
            -fl, A, s_elip

        Returns
        -------
        f2 [deg]
        l2 [deg]
        A2 [deg]

        '''
        f_st = []; l_st = []; A_st =[]
        f_ost = []; l_ost = []; A_ost=[]
        i = 0
        
        if self.f != [''] and self.l != [''] and self.A != [''] and self.s_elip != ['']:
            
            while i < len(self.f):
                
                f_2,l_2,A_2 = self.kivioj(self.f[i], self.l[i], self.A[i], self.s_elip[i], self.a, self.e2)[0:3]
                f_st.append(self.dms(f_2)); l_st.append(self.dms(l_2)); A_st.append(self.dms(A_2))
                f_ost.append(np.rad2deg(f_2)); l_ost.append(np.rad2deg(l_2)); A_ost.append(np.rad2deg(A_2))
                i += 1
                
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [wprost]:' + Format.normal + '\n\u03C62: ',('{} '*len(f_st)).format(*f_st), '\n\u03BB2: ',('{} '*len(l_st)).format(*l_st), '\nA2: ',('{} '*len(A_st)).format(*A_st))
        return(f_ost,l_ost,A_ost)
     
    def odwrotne(self):
        '''
        Wykonuje zadanie odwrotne przy uzyciu ALGORYTMU VINCENTEGO na podstawie dostepnych danych:
            -fl, fl2
            -XYZ, XYZ2

        Returns
        -------
        A [deg]
        A2 [deg]
        s_elip [m]

        '''
        A_st = []; A2_st = []
        A_ost = []; A2_ost = []
        s_elip_ost = []
        
        if self.X != [''] and self.Y != [''] and self.Z != [''] and self.X2 != [''] and self.Y2 != [''] and self.Z2 != ['']:
            
            i = 0
            self.f = []; self.l = []; self.h = []
            self.f2 = []; self.l2 = []; self.h2 = []
            
            while i < len(self.X):
                
                f,l,h = self.xyz2flh(self.X[i], self.Y[i], self.Z[i], self.a, self.e2)
                f2,l2,h2 = self.xyz2flh(self.X2[i], self.Y2[i], self.Z2[i], self.a, self.e2)
                self.f.append(f)
                self.l.append(l)
                self.h.append(h)
                self.f2.append(f2)
                self.l2.append(l2)
                self.h2.append(h2)
                i += 1
                
        if self.f != [''] and self.l != [''] and self.f2 != [''] and self.l2 != ['']:
            
            i = 0
            
            while i < len(self.f):
                
                s_elip, A, A2 = self.vincenty(self.f[i], self.l[i], self.f2[i], self.l2[i], self.a, self.e2)
                A_st.append(self.dms(A)); A2_st.append(self.dms(A2))
                A_ost.append(np.rad2deg(A)); A2_ost.append(np.rad2deg(A2))
                s_elip_ost.append(s_elip)
                i += 1
        
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [odwrotne]:' + Format.normal + '\nA1: ',('{} '*len(A_st)).format(*A_st), '\nA2: ',('{} '*len(A2_st)).format(*A2_st), '\ns_elip: ',('{:.3f} '*len(s_elip_ost)).format(*s_elip_ost), '[m]')
        return(A_ost,A2_ost,s_elip_ost)  
     
    def neu(self):
        '''
        Przelicza wspl pkt 1 na wspl pkt 2 w ukl flh oraz XYZ,
        przy uzyciu ukl wspl topocentrycznych NEU na podstawie dostepnych danych:
            -flh, saz
            -XYZ, saz

        Returns
        -------
        X [m]
        Y [m]
        Z [m]
        f [deg]
        l [deg]
        h [m]
        

        '''
        
        f2_st = []; l2_st = []
        f2_ost = []; l2_ost = []
        h2_ost = []
        X2_ost = []; Y2_ost = []; Z2_ost = []
        i = 0
        
        if self.f != [''] and self.l != [''] and self.h != ['']:
            
            while i < len(self.A):
                
                X,Y,Z = self.flh2xyz(self.f[i], self.l[i], self.h[i], self.a, self.e2)
                DeltaNEU = self.saz2neu(self.s[i], self.A[i], self.z[i])
                DeltaXYZ = self.neu2xyz(DeltaNEU, self.f[i], self.l[i])
                X2 = X + DeltaXYZ[0]
                Y2 = Y + DeltaXYZ[1]
                Z2 = Z + DeltaXYZ[2]
                f2,l2,h2 = self.xyz2flh(X2, Y2, Z2, self.a, self.e2)
                f2_st.append(self.dms(f2)); l2_st.append(self.dms(l2))
                f2_ost.append(np.rad2deg(f2)); l2_ost.append(np.rad2deg(l2))
                h2_ost.append(h2)
                X2_ost.append(X2); Y2_ost.append(Y2); Z2_ost.append(Z2)
                i += 1
                
        if self.X != [''] and self.Y != [''] and self.Z != ['']:
            
            while i < len(self.A):
                
                f,l,h = self.xyz2flh(self.X[i], self.Y[i], self.Z[i], self.a, self.e2)
                DeltaNEU = self.saz2neu(self.s[i], self.A[i], self.z[i])
                DeltaXYZ = self.neu2xyz(DeltaNEU, f, l)
                X2 = self.X[i] + DeltaXYZ[0]
                Y2 = self.Y[i] + DeltaXYZ[1]
                Z2 = self.Z[i] + DeltaXYZ[2]
                f2,l2,h2 = self.xyz2flh(X2, Y2, Z2, self.a, self.e2)
                f2_st.append(self.dms(f2)); l2_st.append(self.dms(l2))
                f2_ost.append(np.rad2deg(f2)); l2_ost.append(np.rad2deg(l2))
                h2_ost.append(h2)
                X2_ost.append(X2); Y2_ost.append(Y2); Z2_ost.append(Z2)
                i += 1

        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [neu]:' + Format.normal + 
              '\n\u03C62: ',('{} '*len(f2_st)).format(*f2_st), '\n\u03BB2: ',('{} '*len(l2_st)).format(*l2_st), '\nh2: ',('{:.3f} '*len(h2_ost)).format(*h2_ost), '[m]'
              '\nX2: ',('{:.3f} '*len(X2_ost)).format(*X2_ost), '[m]\nY2: ',('{:.3f} '*len(Y2_ost)).format(*Y2_ost), '[m]\nZ2: ',('{:.3f} '*len(Z2_ost)).format(*Z2_ost), '[m]')
        return(f2_ost,l2_ost,h2_ost,X2_ost,Y2_ost,Z2_ost)
               
    def Azymut(self):
        '''
        Oblicza wartosc azymutu wprost i odwrotnego dla 
        dostepnych wspolrzednych na plaszczyznie:
           -PL1992
           -PL2000
           -GK, l0
           
        Returns
        -------
        Aab - azymut wprost [deg]
        Aba - azymut odwrotny [deg]
        -------
        alfa ab - azymut wprost na plaszczyznie GK (azymut topograficzny, kat kierunkowy) [deg]
        alfa ba - azymut odwrotny na plaszczyznie GK (azymut topograficzny, kat kierunkowy) [deg]
        gamma a - zbieznosc poludnikow w pkt poczatkowym [deg]
        gamma b - zbieznosc poludnikow w pkt koncowym [deg]
        delta ab - redukcja kierunku wprost [deg]
        delta ba - redukcja kierunku odwrotnego [deg]

        '''
        Aab_ost = []; Aba_ost = []
        Aab_st = []; Aba_st = []
        alfa_ab_st = []; alfa_ba_st = []; gamma_a_st = []; gamma_b_st = []; delta_ab_st = []; delta_ba_st = []
        l0=[]
        self.xgk=[]
        self.ygk=[]
        self.xgk2=[]
        self.ygk2=[]
        i = 0
        
        if self.x1992 != [''] and self.y1992 != [''] and self.x1992_2 != [''] and self.y1992_2 != ['']:
            
            while i < len(self.x1992):
                #do zmiany/popr
                self.xgk.append((self.x1992[i] + 5300000)/self.m0_1992)
                self.ygk.append((self.y1992[i] - 500000)/self.m0_1992)                
                self.xgk2.append((self.x1992_2[i] + 5300000)/self.m0_1992)
                self.ygk2.append((self.y1992_2[i] - 500000)/self.m0_1992)
                l0.append(np.radians(19))
                i += 1
                
        elif self.x2000 != [''] and self.y2000 != [''] and self.x2000_2 != [''] and self.y2000_2 != ['']:
            
            while i < len(self.x2000):
                #do zmiany/popr
                ns1 = self.strefa2(self.y2000[i])
                ns2 = self.strefa2(self.y2000_2[i])
                self.xgk.append(self.x2000[i]/self.m0_2000)
                self.ygk.append((self.y2000[i] - 500000 -  ns1*1000000)/self.m0_2000)
                self.xgk2.append(self.x2000_2[i]/self.m0_2000)
                self.ygk2.append((self.y2000_2[i] - 500000 - ns2*1000000)/self.m0_2000)
                if ns1 == ns2:
                    if ns1 == 5:
                        l0.append(np.radians(15))
                    elif ns1 == 6:
                        l0.append(np.radians(18))
                    elif ns1 == 7:
                        l0.append(np.radians(21))
                    elif ns1 == 8:
                        l0.append(np.radians(24))   
                i += 1
               
        if self.xgk != [] and self.ygk != [] and self.xgk2 != [] and self.ygk2 != []:
            
            i = 0
            
            while i < len(self.xgk):

                #krok 2 - policzenie azymutu (kierunku) na plaszczyznie gk:
                alfa_ab = np.arctan2(self.ygk2[i]-self.ygk[i],self.xgk2[i]-self.xgk[i])
                alfa_ba = np.arctan2(self.ygk[i]-self.ygk2[i],self.xgk[i]-self.xgk2[i])
                alfa_ab_st.append(self.dms(alfa_ab)); alfa_ba_st.append(self.dms(alfa_ba))
                #krok 3 - policzenie zbieznosci pld:
                gamma_a = self.zbiez_pld(self.xgk[i], self.ygk[i], self.a, self.e2)
                gamma_b = self.zbiez_pld(self.xgk2[i], self.ygk2[i], self.a, self.e2)
                gamma_a_st.append(self.dms(gamma_a)); gamma_b_st.append(self.dms(gamma_b))
                #krok 4 - policzenie redukcji kierunku:
                # l0 ---> idk co z tym zrobic jak bd podane wspl w gk
                delta_ab = self.zbiez_kier(self.xgk[i], self.ygk[i], self.xgk2[i], self.ygk2[i], self.a, self.e2, l0[i])
                delta_ba = self.zbiez_kier(self.xgk2[i], self.ygk2[i], self.xgk[i], self.ygk[i], self.a, self.e2, l0[i])
                delta_ab_st.append(self.dms(delta_ab)); delta_ba_st.append(self.dms(delta_ba))
                
                A_ab = alfa_ab + delta_ab + gamma_a
                A_ba = alfa_ba + delta_ba + gamma_b
                
                if A_ab < 0:
                    A_ab = A_ab + 2*np.pi
                if A_ba < 0:
                    A_ba = A_ba + 2*np.pi 
                    
                Aab_st.append(self.dms(A_ab)), Aba_st.append(self.dms(A_ba))
                Aab_ost.append(np.rad2deg(A_ab)), Aba_ost.append(np.rad2deg(A_ba))

                i += 1
        
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [Azymut]:' + Format.normal)     
        if self.posrednie == True:
            print(Format.kursywa + 'Wyniki posrednie:' + Format.normal + '\n\u03B1ab: ',('{} '*len(alfa_ab_st)).format(*alfa_ab_st) + '\n\u03B1ba: ',('{} '*len(alfa_ba_st)).format(*alfa_ba_st) + '\n\u03B3a: ',('{} '*len(gamma_a_st)).format(*gamma_a_st)  + '\n\u03B3b: ',('{} '*len(gamma_b_st)).format(*gamma_b_st) + '\n\u03B4ab: ',('{} '*len(delta_ab_st)).format(*delta_ab_st) + '\n\u03B4ba: ',('{} '*len(delta_ba_st)).format(*delta_ba_st) + Format.kursywa + '\nWyniki ostateczne:' + Format.normal)
        print('Aab: ',('{} '*len(Aab_st)).format(*Aab_st),'\nAba: ',('{} '*len(Aba_st)).format(*Aba_st))
        return(Aab_ost,Aba_ost)
    
    def OdlElip(self):
        '''
        Oblicza odleglosc na elipsoidzie (dlugosc lini geodezyjnej) dla 
        dostepnych wspolrzednych na plaszczyznie:
            -PL1992
            -PL2000
            -GK, l0
            
        
        Returns
        -------
        s_elip - dlugosc na elipsoidzie [m]
        -------
        s_gk - dlugosc na plszczyznie GK [m]
        r_gk - redukcja odleglosci na plaszczyzne GK [m]

        '''
        s_elip_ost = []; s_gk_ost=[]; r_gk_ost = []
        l0 =[]
        self.xgk=[]
        self.ygk=[]
        self.xgk2=[]
        self.ygk2=[]
        i = 0
        
        if self.x1992 != [''] and self.y1992 != [] and self.x1992_2 != [''] and self.y1992_2 != ['']:
            
            while i < len(self.x1992):
                #do zmiany/popr
                self.xgk.append((self.x1992[i] + 5300000)/self.m0_1992)
                self.ygk.append((self.y1992[i] - 500000)/self.m0_1992)                
                self.xgk2.append((self.x1992_2[i] + 5300000)/self.m0_1992)
                self.ygk2.append((self.y1992_2[i] - 500000)/self.m0_1992)
                l0.append(np.radians(19))
                i += 1
                
        elif self.x2000 != [''] and self.y2000 != [] and self.x2000_2 != [''] and self.y2000_2 != ['']:
            
            
            while i < len(self.x2000):
                #do zmiany/popr
                ns1 = self.strefa2(self.y2000[i])
                ns2 = self.strefa2(self.y2000_2[i])
                self.xgk.append(self.x2000[i]/self.m0_2000)
                self.ygk.append((self.y2000[i] - 500000 -  ns1*1000000)/self.m0_2000)
                self.xgk2.append(self.x2000_2[i]/self.m0_2000)
                self.ygk2.append((self.y2000_2[i] - 500000 - ns2*1000000)/self.m0_2000)
                if ns1 == ns2:
                    if ns1 == 5:
                        l0.append(np.radians(15))
                    elif ns1 == 6:
                        l0.append(np.radians(18))
                    elif ns1 == 7:
                        l0.append(np.radians(21))
                    elif ns1 == 8:
                        l0.append(np.radians(24))   
                i += 1
        if self.xgk != [] and self.ygk != [] and self.xgk2 != [] and self.ygk2 != []:  
            
            i = 0
            
            while i < len(self.xgk):
                #krok 1: s_gk
                # s_gk = np.sqrt((self.xgk2[i] - self.xgk[i])**2 + (self.ygk2[i] - self.ygk[i])**2)
                #krok 2: r_ab
                r, s_gk = self.red_gk(self.xgk[i], self.ygk[i], self.xgk2[i], self.ygk2[i], l0[i], self.a, self.e2)
                #krok 3: s_elip
                s_elip = s_gk - r
                s_elip_ost.append(s_elip); s_gk_ost.append(s_gk); r_gk_ost.append(r)
                
                i += 1
                
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [OdlElip]:' + Format.normal)     
        if self.posrednie == True:
            print(Format.kursywa + 'Wyniki posrednie:' + Format.normal + '\ns_gk: ',('{:.3f} '*len(s_gk_ost)).format(*s_gk_ost) + '[m]\nr_gk: ',('{:.3f} '*len(r_gk_ost)).format(*r_gk_ost) + '[m]' + Format.kursywa + '\nWyniki ostateczne:' + Format.normal)
        print('s_elip: ',('{:.3f} '*len(s_elip_ost)).format(*s_elip_ost) + '[m]')
        return(s_elip_ost)

    def OdlPL2000(self):
        '''
        Oblicza odleglosc w układzie PL-2000 dla 
        dostepnych danych:
            -PL2000
            -XYZ, XYZ2
            -fl, fl2

        Returns
        -------
        s_2000 - dlugosc w ukladzie PL2000 [m]
        -------
        s_XYZ - dlugosc w ukladzie wspl prostokatnych [m]
        s_0 - dlugosc skosna [m]
        s_elip - dlugosc na elipsoidzie [m]
        s_gk - dlugosc na plszczyznie GK [m]
        r_gk - redukcja odleglosci na plaszczyzne GK [m]
        
        '''
        l0 = []
        s_2000_ost = []; s_elip_ost = []; s_0_ost = []; s_gk_ost = []; r_gk_ost = []; s_XYZ_ost = []
        i = 0
        
        if self.x2000 != [''] and self.y2000 != [''] and self.x2000_2 != [''] and self.y2000_2 != ['']:

            while i < len(self.x2000):
                
                s_2000 = np.sqrt((self.x2000_2[i] - self.x2000[i])**2 + (self.y2000_2[i] - self.y2000[i])**2)
                s_2000_ost.append(s_2000)
                i += 1
        
        else:
            
            if self.X != [''] and self.Y != [''] and self.Z != [''] and self.X2 != [''] and self.Y2 != [''] and self.Z2 != ['']:
                
                self.f = []
                self.l = []
                self.h = []
                self.f2 = []
                self.l2 = []
                self.h2 = []
                i = 0
                
                while i < len(self.X):
                    
                    f,l,h = self.xyz2flh(self.X[i], self.Y[i], self.Z[i], self.a, self.e2)
                    self.f.append(f)
                    self.l.append(l)
                    self.h.append(h)
                    f2,l2,h2 = self.xyz2flh(self.X2[i], self.Y2[i], self.Z2[i], self.a, self.e2)
                    self.f2.append(f2)
                    self.l2.append(l2)
                    self.h2.append(h2)
                    
                    i += 1
                    
            if self.f != [''] and self.l != [''] and self.f2 != [''] and self.l2 != ['']:
                
                    i = 0  
                    
                    while i < len(self.f):     
                       
                        ns = self.strefa(self.l[i])
                        x2000,y2000,xgk,ygk = self.fl2PL2000(self.f[i], self.l[i], self.a, self.e2, ns)
                        ns2 = self.strefa(self.l2[i])
                        x2000_2,y2000_2,xgk_2,ygk_2 = self.fl2PL2000(self.f2[i], self.l2[i], self.a, self.e2, ns2)
                        
                        if ns == ns2:
                            if ns == 5:
                                l0.append(np.radians(15))
                            elif ns == 6:
                                l0.append(np.radians(18))
                            elif ns == 7:
                                l0.append(np.radians(21))
                            elif ns == 8:
                                l0.append(np.radians(24)) 
                        
                        s_XYZ = np.sqrt((self.X2[i]-self.X[i])**2 + (self.Y2[i]-self.Y[i])**2 + (self.Z2[i]-self.Z[i])**2)
                        s_0, s_elip = self.redu_d(xgk, xgk_2, ygk, ygk_2, h, h2, s_XYZ, l0[i], self.a, self.e2) 
                        r_gk, s_gk = self.red_gk(xgk, ygk, xgk_2, ygk_2, l0[i], self.a, self.e2)
                        s_gk = s_elip + r_gk
                        s_2000 = 0.999923*s_gk
                        
                        s_2000_ost.append(s_2000); s_elip_ost.append(s_elip); s_0_ost.append(s_0); s_gk_ost.append(s_gk); r_gk_ost.append(r_gk); s_XYZ_ost.append(s_XYZ)
                        i += 1
                
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [OdlPL2000]:' + Format.normal)     
        if self.posrednie == True:
            print(Format.kursywa + 'Wyniki posrednie:' + Format.normal + '\ns_XYZ: ',('{:.3f} '*len(s_XYZ_ost)).format(*s_XYZ_ost) + '[m]\ns_0: ',('{:.3f} '*len(s_0_ost)).format(*s_0_ost) + '[m]\ns_elip: ',('{:.3f} '*len(s_elip_ost)).format(*s_elip_ost) + '[m]\ns_gk: ',('{:.3f} '*len(s_gk_ost)).format(*s_gk_ost) + '[m]\nr_gk: ',('{:.3f} '*len(r_gk_ost)).format(*r_gk_ost) + '[m]' + Format.kursywa + '\nWyniki ostateczne:' + Format.normal)
        print('s_2000: ',('{:.3f} '*len(s_2000_ost)).format(*s_2000_ost) + '[m]')
        return(s_2000)
    
    def OdlPL1992(self):
        pass
    
    def OdlSkosna(self):
        pass
    
    def OdlNEU(self):
        pass
    
    def Znieksztalcenia(self):
        '''
        Oblicza znieksztalcenia pola powierzchni i odleglosci dla wspl w podanym
        ukl wspl plaskich:
            -PL2000
            -PL1992

        Returns
        -------
        znie_odl - znieksztalcenie odlegosci [cm/km]
        znie_pol - znieksztalcenie pola powierzchni [m2/ha]
        -------
        mgk - skala ukl GK [-]
        m2000 - skala ukl PL2000 [-]
        m1992 - skala ukl PL1992 [-]

        '''
        znie_odl_ost = []; znie_pol_ost = []
        mgk_ost = []; m2000_ost = []; m1992_ost =[]
        if self.x2000 != [''] and self.y2000 != ['']:
            
            i = 0
            
            while i < len(self.x2000): 
                
                ns = self.strefa2(self.y2000[i])
                f,l,xgk,ygk = self.PL20002fl(self.x2000[i], self.y2000[i], self.a, self.e2, ns)
                R = np.sqrt(self.Np(f, self.a, self.e2)*self.Mp(f, self.a, self.e2))
                mgk = 1 + ygk**2/2/R**2 + ygk**4/24/R**4
                m2000 = self.m0_2000 * mgk
                znie_odl = 1000 * m2000 - 1000
                znie_pol = 10000 * m2000**2 - 10000
                znie_odl_ost.append(znie_odl); znie_pol_ost.append(znie_pol)
                mgk_ost.append(mgk); m2000_ost.append(m2000)
                i += 1
                
        elif self.x1992 != [''] and self.y1992 != ['']:
             
             i = 0
             
             while i < len(self.x1992): 
                 
                 f,l,xgk,ygk = self.PL19922fl(self.x1992[i], self.y1992[i], self.a, self.e2)
                 R = np.sqrt(self.Np(f, self.a, self.e2)*self.Mp(f, self.a, self.e2))
                 mgk = 1 + ygk**2/2/R**2 + ygk**4/24/R**4
                 m1992 = self.m0_1992 * mgk
                 znie_odl = 1000 * m1992 - 1000
                 znie_pol = 10000 * m1992**2 - 10000
                 znie_odl_ost.append(znie_odl); znie_pol_ost.append(znie_pol)
                 mgk_ost.append(mgk), m1992_ost.append(m1992)
                 i += 1     
                 
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [Znieksztalcenia]:' + Format.normal)     
        if self.posrednie == True:
            print(Format.kursywa + 'Wyniki posrednie:' + Format.normal + '\nmgk: ',('{:.5f} '*len(mgk_ost)).format(*mgk_ost) + '[-]\nm2000: ',('{:.3f} '*len(m2000_ost)).format(*m2000_ost) + '[-]\nm1992: ',('{:.3f} '*len(m1992_ost)).format(*m1992_ost) + '[-]' + Format.kursywa + '\nWyniki ostateczne:' + Format.normal)
        print('znieksztalcenie odleglosci: ',('{:.3f} '*len(znie_odl_ost)).format(*znie_odl_ost) + '[cm/km]\nznieksztalcenie pola: ',('{:.3f} '*len(znie_pol_ost)).format(*znie_pol_ost) + '[m2/ha]')
        return(znie_odl_ost,znie_pol_ost)

    def Transformacja(self):
        '''
        Wykonuje transformacje współrzędnych pierwotnych XYZ do układu wtótnego XYZ2 
        na podstawie dostępnych danych. Model transformacji wybierany jest 
        na podstawie ilości dostarczonych parametrów transformacji p:
            -XYZ, p = [kappa x, kappa y, kappa z, alfa, beta, gamma, dx, dy, dz] 
                -> transformacja kwaziafiniczna 
            -XYZ, p = [kappa, alfa, beta, gamma, dx, dy, dz] 
                -> transformacja Bursa-Wolfa
            -XYZ, p = [alfa, beta, gamma, dx, dy, dz]
                -> transformacja izometryczna

        Returns
        -------
            -XYZ2 - wspl prostokatne w ukladzie wtornym [m]

        '''
        
        X2_ost = []; Y2_ost = []; Z2_ost = []
        i = 0
        
        if self.X != [''] and self.Y != [''] and self.Z != ['']:
            
            while i < len(self.X):

                if len(self.p[i]) == 9:

                    X2, Y2, Z2 = self.kwazi([self.X[i], self.Y[i], self.Z[i]], self.p[i])
                    X2_ost.append(X2), Y2_ost.append(Y2), Z2_ost.append(Z2)
                    
                if len(self.p[i]) == 7:
                    
                    X2, Y2, Z2 = self.bursa([self.X[i], self.Y[i], self.Z[i]], self.p[i])
                    X2_ost.append(X2), Y2_ost.append(Y2), Z2_ost.append(Z2)
                    
                if len(self.p[i]) == 6:
                    
                    X2, Y2, Z2 = self.izometr([self.X[i], self.Y[i], self.Z[i]], self.p[i])
                    X2_ost.append(X2), Y2_ost.append(Y2), Z2_ost.append(Z2)
                
                i += 1
                
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [Transformacje]:' + Format.normal)     
        if self.posrednie == True:
            pass
            #print(Format.kursywa + 'Wyniki posrednie:' + Format.normal + '\ns_XYZ: ',('{:.3f} '*len(s_XYZ_ost)).format(*s_XYZ_ost) + '[m]\ns_0: ',('{:.3f} '*len(s_0_ost)).format(*s_0_ost) + '[m]\ns_elip: ',('{:.3f} '*len(s_elip_ost)).format(*s_elip_ost) + '[m]\ns_gk: ',('{:.3f} '*len(s_gk_ost)).format(*s_gk_ost) + '[m]\nr_gk: ',('{:.3f} '*len(r_gk_ost)).format(*r_gk_ost) + '[m]' + Format.kursywa + '\nWyniki ostateczne:' + Format.normal)                        
        print('X2: ',('{:.3f} '*len(X2_ost)).format(*X2_ost), '[m]\nY2: ',('{:.3f} '*len(Y2_ost)).format(*Y2_ost), '[m]\nZ2: ',('{:.3f} '*len(Z2_ost)).format(*Z2_ost), '[m]')
        return(X2_ost,Y2_ost,Z2_ost) 
      
if __name__=='__main__':
      
    # T E S T Y 
    
    proba1 = skrypt(x1992=100, y1992=7400000)
    proba1.flh()

    proba2 = skrypt(x2000=100, y2000=7400000, h=100)
    proba2.XYZ()
    
    proba3 = skrypt(x1992=463717.558, y1992=294552.995, h=100)
    proba3.PL2000()
    
    proba4 = skrypt(s_elip=43000.0,A=230,f='54 7 20.79937',l='23 0 26.12508')
    proba4.wprost()
    
    proba5 = skrypt(f='54 7 20.79937',l='23 0 26.12508',f2='53 52 23.05857',l2='22 30 23.24978')
    proba5.odwrotne()
    
    proba6 = skrypt(x2000=5763554.505, y2000=5569082.651)
    proba6.PL1992()
    
    proba7 = skrypt(x2000 = 5906044.038,y2000 = 6481105.250,x2000_2 = 5895976.242,y2000_2 = 6498384.639)
    proba7.Azymut()
    
    proba8 = skrypt(p=[0.000002, 0.000001, -0.000002, '-0 0 2', '0 0 0.5', '0 0 1', 100, 200, -300],X=5500000, Y=4350000,Z=3125000)
    proba8.Transformacja()

    