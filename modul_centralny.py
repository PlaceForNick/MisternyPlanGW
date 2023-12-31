import numpy as np
from math import *
from modul_obliczeniowy import *

class Format:
    normal = '\033[0m'
    podkresl = '\033[4m'
        
class skrypt(funkcje):
    
    ID = 0
    
    def __init__(self, model='grs80', zapis=False, posrednie=False, nazwa='', X='', Y='', Z='', f='', l='', h='', X2='', Y2='', Z2='', s='', A='', z ='', x2000='', y2000='', x1992='', y1992='', xgk='', ygk='', f2='', l2='', h2='', s_elip='', x2000_2='', y2000_2='', x1992_2='', y1992_2='', xgk2='', ygk2=''):

        skrypt.ID += 1
        
        self.__elipsoida(model) #wybor elipsoidy
        # self.__zapiszplik(zapis, nazwa) #wybor zapisu do pliku txt     
        self.posrednie = posrednie #definiuje czy ma wypluwac obliczenia posrednie
        
        self.m0_1992 = 0.9996
        self.m0_2000 = 0.999923
        
        #zamiana podanych danych na liste
        dane = [X, Y, Z, f, l, h, X2, Y2, Z2, s, A, z, x2000, y2000, x1992, y1992, xgk, ygk, f2, l2, h2, s_elip, x2000_2, y2000_2, x1992_2, y1992_2, xgk2, ygk2]
        dane_ost = []
        for wartosc in dane:
            if type(wartosc) == list:
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
        
        #zamiana stopni na radiany           
        dane_kat = [self.f, self.l, self.A, self.z, self.f2, self.l2]
        dane_kat_ost = []
        for wartosc_lista in dane_kat:
            wartosc_ost = []
            i = 0  
            while True:
                try:
                    wartosc = wartosc_lista[i]
                    if type(wartosc) == list:
                        wartosc = wartosc[0]
                    if wartosc =='':
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
                
                f,l,xgk,ygk = self.PL20002fl(self.x2000[i], self.y2000[i], self.a, self.e2, self.strefa2(self.y2000[i])[0])
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

        Returns
        -------
        A [deg]
        A2 [deg]
        s_elip [m]

        '''
        A_st = []; A2_st = []
        A_ost = []; A2_ost = []
        s_elip_ost = []
        i = 0
        
        if self.f != [''] and self.l != [''] and self.f2 != [''] and self.l2 != ['']:
            
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
        Aab_ost = []; Aba_ost = []
        Aab_st = []; Aba_st = []
        l0=[]
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
                l0.append(''),self.xgk.append(''),self.xgk.append(''),self.xgk2.append(''),self.ygk2.append('')
                ns, l0[i] = self.strefa2(self.y2000[i])
                self.xgk[i], self.ygk[i] = (self.PL20002gk(self.x2000[i], self.y2000[i], ns))
                ns, l0[i] = self.strefa2(self.y2000_2[i])
                self.xgk2[i], self.ygk2[i] = (self.PL20002gk(self.x2000_2[i], self.y2000_2[i], ns))
                
                i += 1
               
        if self.xgk != [] and self.ygk != [] and self.xgk2 != [] and self.ygk2 != []:
        # if True:
            
            i = 0
            
            while i < len(self.xgk):

                #krok 2 - policzenie azymutu (kierunku) na plaszczyznie gk:
                alfa_ab = np.arctan2(self.ygk2[i]-self.ygk[i],self.xgk2[i]-self.xgk[i])
                alfa_ba = np.arctan2(self.ygk[i]-self.ygk2[i],self.xgk[i]-self.xgk2[i])
                # print(self.dms(alfa_ab),self.dms(alfa_ba))
                #krok 3 - policzenie zbieznosci pld:
                gamma_a = self.zbiez_pld(self.xgk[i], self.ygk[i], self.a, self.e2)
                gamma_b = self.zbiez_pld(self.xgk2[i], self.ygk2[i], self.a, self.e2)
                # print(self.dms(gamma_a), self.dms(gamma_b))
                #krok 4 - policzenie redukcji kierunku:
                # l0 ---> idk co z tym zrobic jak bd podane wspl w gk
                delta_ab = self.zbiez_kier(self.xgk[i], self.ygk[i], self.xgk2[i], self.ygk2[i], self.a, self.e2, l0[i])
                delta_ba = self.zbiez_kier(self.xgk2[i], self.ygk2[i], self.xgk[i], self.ygk[i], self.a, self.e2, l0[i])
                # print(self.dms(delta_ab), self.dms(delta_ba))
                A_ab = alfa_ab + delta_ab + gamma_a
                A_ba = alfa_ba + delta_ba + gamma_b
                
                if A_ab < 0:
                    A_ab = A_ab + 2*np.pi
                if A_ba < 0:
                    A_ba = A_ba + 2*np.pi 
                    null
                Aab_st.append(self.dms(A_ab)), Aba_st.append(self.dms(A_ba))
                Aab_ost.append(np.rad2deg(A_ab)), Aba_ost.append(np.rad2deg(A_ba))

                i += 1
                
        print(Format.podkresl + f'\nZapytanie {skrypt.ID} [Azymut]:' + Format.normal + '\nAab: ',('{} '*len(Aab_st)).format(*Aab_st),'\nAba: ',('{} '*len(Aba_st)).format(*Aba_st))# '\nA2: ',('{} '*len(A2_st)).format(*A2_st), '\ns_elip: ',('{:.3f} '*len(s_elip_ost)).format(*s_elip_ost), '[m]')
        return(Aab_ost,Aba_ost)
    
    def odl_elip(self):
        
        s_gk = np.sqrt((self.xgk2[i] - self.xgk[i])**2 + (self.ygk2[i]-self.ygk[i])**2)
        
        # red_ab = 
        
        s_elip = s_gk - r
    
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

    