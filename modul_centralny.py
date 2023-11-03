import numpy as np
from math import *
from modul_obliczeniowy import *

class skrypt(funkcje):
    
    def __init__(self, model='grs80', zapis=False, posrednie=False, nazwa='', X='', Y='', Z='', f='', l='', h='', X2='', Y2='', Z2='', s_elip='', A='', z ='', x2000='', y2000='', x1992='', y1992='', xgk='', ygk=''):
        
        self.__elipsoida(model) #wybor elipsoidy
        # self.__zapiszplik(zapis, nazwa) #wybor zapisu do pliku txt     
        self.posrednie = posrednie #definiuje czy ma wypluwac obliczenia posrednie
        
        #zamiana podanych danych na liste
        dane = [X, Y, Z, f, l, h, X2, Y2, Z2, s_elip, A, z, x2000, y2000, x1992, y1992, xgk, ygk]
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
        self.s_elip = dane_ost[9]
        self.A = dane_ost[10]
        self.z = dane_ost[11]
        self.x2000 = dane_ost[12]
        self.y2000 = dane_ost[13]
        self.x1992 = dane_ost[14]
        self.y1992 = dane_ost[15]
        self.xgk = dane_ost[16]
        self.ygk = dane_ost[17]
        
        #zamiana stopni na radiany           
        dane_kat = [self.f, self.l, self.A, self.z]
        dane_kat_ost = []
        for wartosc_lista in dane_kat:
            wartosc_ost = []
            i = 0  
            while True:
                try:
                    wartosc = wartosc_lista[i]
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
         
        f_st = []; l_st = []
        h_ost = []
        i = 0   
        
        if self.X != [''] and self.Y != [''] and self.Z != ['']:
            
            while i < len(self.X):
                
                f,l,h = self.xyz2flh(self.X[i], self.Y[i], self.Z[i], self.a, self.e2)
                f_st.append(self.dms(f)); l_st.append(self.dms(l))
                h_ost.append(h)
                i += 1
        
        elif self.x2000 != [''] and self.y2000 != ['']:
            
            while i < len(self.x2000):

                ns = int(str(self.y2000[i])[0])
                f,l,xgk,ygk = self.PL20002fl(self.x2000[i], self.y2000[i], self.a, self.e2, ns)
                f_st.append(self.dms(f)); l_st.append(self.dms(l))
                h_ost = self.h
                i += 1
        
        elif self.x1992 != [''] and self.y1992 != ['']:
            
            while i < len(self.x1992):

                f,l,xgk,ygk = self.PL19922fl(self.x1992[i], self.y1992[i], self.a, self.e2)
                f_st.append(self.dms(f)); l_st.append(self.dms(l))
                h_ost = self.h
                i += 1

        print('\n\u03C6: ',('{} '*len(f_st)).format(*f_st), '\n\u03BB: ',('{} '*len(l_st)).format(*l_st), '\nh: ',('{} '*len(h_ost)).format(*h_ost), '[m]')
        
    def XYZ(self):
        
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
                
        print( '\nX: ',('{:.3f} '*len(X_ost)).format(*X_ost), '[m]\nY: ',('{:.3f} '*len(Y_ost)).format(*Y_ost), '[m]\nZ: ',('{:.3f} '*len(Z_ost)).format(*Z_ost), '[m]' )

    def PL2000(self):
        
        x_ost = []; y_ost = []
        i = 0
        
        if self.f != [''] and self.l != ['']:
            
            while i < len(self.f):
                
                if self.l[i] < radians(16.5) and self.l[i] > radians(13.5): #ns = 5
                    l0 = radians(15)
                    ns = 5
                elif self.l[i] < radians(19.5) and self.l[i] > radians(16.5): #ns = 6
                    l0 = radians(18)
                    ns = 6
                elif self.l[i] < radians(22.5) and self.l[i] > radians(19.5): #ns = 7
                    l0 = radians(21)
                    ns = 7
                elif self.l[i] < radians(25.5) and self.l[i] > radians(22.5): #ns = 8
                    l0 = radians(24)
                    ns = 8
                    
                x2000,y2000,xgk,ygk = self.fl2PL2000(self.f[i], self.l[i], self.a, self.e2, ns)
                x_ost.append(x2000); y_ost.append(y2000)
                i += 1
                
        elif self.X != [''] and self.Y != [''] and self.Z != ['']:
            
            while i < len(self.X):
                
                f,l,h = self.xyz2flh(self.X[i], self.Y[i], self.Z[i], self.a, self.e2)
                
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
                 
                x2000,y2000,xgk,ygk = self.fl2PL2000(f, l, self.a, self.e2, ns)
                x_ost.append(x2000); y_ost.append(y2000)
                i += 1
                
        elif self.x1992 != [''] and self.y1992 != ['']:
            
            while i < len(self.x1992):
                
                f,l,xgk,ygk = self.PL19922fl(self.x1992[i], self.y1992[i], self.a, self.e2)
                
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
                 
                x2000,y2000,xgk,ygk = self.fl2PL2000(f, l, self.a, self.e2, ns)
                x_ost.append(x2000); y_ost.append(y2000)
                i += 1
                
        print( '\nx2000: ',('{:.3f} '*len(x_ost)).format(*x_ost), '[m]\ny2000: ',('{:.3f} '*len(y_ost)).format(*y_ost), '[m]' )
        
    def wprost(self): #algorytm kivioja
        
        f_st = []; l_st = []; A_st =[]
        i = 0
        
        if self.f != [''] and self.l != [''] and self.A != [''] and self.s_elip != ['']:
            
            while i < len(self.f):
                
                f_2,l_2,A_2 = self.kivioj(self.f[i], self.l[i], self.A[i], self.s_elip[i], self.a, self.e2)[0:3]
                f_st.append(self.dms(f_2)); l_st.append(self.dms(l_2)); A_st.append(self.dms(A_2))
                i += 1
                
        print( '\n\u03C62: ',('{} '*len(f_st)).format(*f_st), '\n\u03BB2: ',('{} '*len(l_st)).format(*l_st), '\nA2: ',('{} '*len(A_st)).format(*A_st))

     
    def odwrotne(self): #algorytm vincentego      
        
        A_st = []; A_2_st = []
        s_elip_ost = []
        i = 0
        
        if self.f != [''] and self.l != [''] and self.f_2 != [''] and self.l_2 != ['']:
            
            while i < len(self.f):
                
                s_elip, A, A_2 = self.vincenty(self.f[i], self.l[i], self.f_2[i], self.l_2[i], self.a, self.e2)
                A_st.append(self.dms(A)); A_2_st.append(self.dms(A_2))
                s_elip_ost.append(s_elip)
                i += 1
        
        print( '\nA1: ',('{} '*len(A_st)).format(*A_st), '\nA2: ',('{} '*len(A_2_st)).format(*A_2_st), '\ns_elip: ',('{:.3f} '*len(s_elip_ost)).format(*s_elip_ost), '[m]')
            
if __name__=='__main__':
    
    proba1 = skrypt(x1992=100, y1992=7400000)
    proba1.flh()
    
    proba2 = skrypt(x2000=100, y2000=7400000, h=100)
    proba2.XYZ()
    
    proba3 = skrypt(x1992=463717.558, y1992=294552.995, h=100)
    proba3.PL2000()
    
    proba4 = skrypt(s_elip=43000.0,A=230,f='54 7 20.79937',l='23 0 26.12508')
    proba4.wprost()