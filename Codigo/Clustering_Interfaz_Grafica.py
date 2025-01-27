#!/usr/bin/env python
# coding: utf-8

# In[24]:


import pandas
import sys
import numpy as np
import random as rd
import math
import matplotlib.pyplot as plt #Imports para dibujar un circulo


####################################################################################################################
############################################### Estructuras de Datos ###############################################
####################################################################################################################


# Parámetros: x,y del punto
class Punto(): 
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
    def print(self):
        print('(',self.x, ', ',self.y,')')
    def __repr__(self):
        return '(%s, %s)' % (self.x,self.y)
    def __str__(self):
        return '(%s, %s)' % (self.x,self.y)
    
    def __eq__(self, other):
        return ((self.x, self.y) == (other.x, other.y))

# Parámetros: x,y (centro), radio
class Circunferencia(): 
    def __init__(self, x, y, radio):
        self.centro = Punto(x,y)
        self.radio = radio
        
    def print(self):
        print('Centro: ',self.centro, ' Radio: ',self.radio)
    def __repr__(self):
        return 'Centro: %s Radio: %s' % (self.centro,self.radio)
    def __str__(self):
        return 'Centro: %s Radio: %s' % (self.centro,self.radio)
    
# Parámetros: circunferencia
class Cluster(): 
    def __init__(self,circunferencia):
        self.circunferencia = circunferencia
        self.puntos = []
        
    def anyadir_punto(self,punto):
        self.puntos.append(punto)
        
    def anyadir_puntos(self,puntos):
        self.puntos = puntos
        
        
    def print(self):
        print('Circunferencia: ',self.circunferencia, ' Puntos: ',self.puntos)
    def __repr__(self):
        return 'Circunferencia: %s Puntos: %s' % (self.circunferencia,self.puntos)
    def __str__(self):
        return 'Circunferencia: %s Puntos: %s' % (self.circunferencia,self.puntos)
    
    

##########################################################################################################################################
############################################### Funciones para leer, crear y mostrar datos ###############################################
##########################################################################################################################################


#Parámetros: Archivo CSV donde se encuentran los puntos
def leer_puntos(archivo_csv):
   
    puntoscsv = pandas.read_csv(archivo_csv, header=None, names=['x', 'y']).values
    puntos = []
    for i in range(len(puntoscsv)):
        puntos.append(Punto(puntoscsv[i][0],puntoscsv[i][1]))
        
    return puntos


# Parámetros: Array de puntos
def dibujar_puntos(puntos):  # learningaboutelectronics.com/Articles/How-to-draw-a-circle-using-matplotlib-in-Python.php
    
    Xs = []
    Ys = []
    
    for i in range(len(puntos)): 
        Xs.append(puntos[i].x)
        Ys.append(puntos[i].y)
        
    plt.plot(Xs, Ys, 'ro', color='black')
    plt.axis('scaled')
    
    plt.savefig("Data/puntosIniciales.png")
    plt.close()
    
    
# Parámetros: Circunferencia con Puntos
def dibujar_clusters(clusters):  # learningaboutelectronics.com/Articles/How-to-draw-a-circle-using-matplotlib-in-Python.php
    
    # Lista de colores con los que se representará las circunferencias con sus puntos asociados
    colores = ['red','blue','green','black','purple','cyan','pink','magenta','orange','brown','yellow']
    
    for i in range(len(clusters)):
       
        cluster = clusters[i]
        Xs = []
        Ys = []
        
        # Dibujamos los puntos asociados a la circunferencia
        for p in cluster.puntos:
            Xs.append(p.x)
            Ys.append(p.y)
        
        plt.plot(Xs, Ys, 'ro', color=colores[i])
        
        #Dibujamos la circunferencia
        circulo = plt.Circle((cluster.circunferencia.centro.x, cluster.circunferencia.centro.y), cluster.circunferencia.radio, color=colores[i%11], fill=False)
        ax=plt.gca()
        ax.add_patch(circulo)
        plt.axis('scaled')
    
    plt.savefig("Data/resultado.png")
    plt.close()


# Parámetros: array de las circunferencias con puntos, redondear = 0: no redondear, 1: redondear, mostrar_puntos = 0: no mostrar, 1: mostrar
def mostrar_clusters(clusters, redondear, mostrar_puntos):
    
    for i in range(len(clusters)):
        
        centrox = clusters[i].circunferencia.centro.x
        centroy = clusters[i].circunferencia.centro.y
        radio = clusters[i].circunferencia.radio
        puntos = clusters[i].puntos
        
        if(redondear):
            centrox = round(centrox, 2)
            centroy = round(centroy, 2)
            radio = round(radio, 2)
        
        print("Circunferencia Nº %d \n" % (i+1))
        print("Centro: (%f, %f), Radio: %f" % (centrox, centroy, radio))
        if(mostrar_puntos):
            print("Puntos: %s \n\n" % puntos)
        else:
            print("\n")
        
        
        
# Parámetros: nombre del archivo resultante, array de las circunferencias con puntos, redondear = 0: no redondear, 1: redondear
def escribir_fichero_clusters(nombre_archivo, clusters, redondear):
     
    f = open("%s.txt" % nombre_archivo,"w+")
        
    for i in range(len(clusters)):
        
        centrox = clusters[i].circunferencia.centro.x
        centroy = clusters[i].circunferencia.centro.y
        radio = clusters[i].circunferencia.radio
        puntos = clusters[i].puntos
        
        if(redondear):
            centrox = round(centrox, 2)
            centroy = round(centroy, 2)
            radio = round(radio, 2)
        
        f.write("Circunferencia %d\n\n" % (i+1))
        f.write("%f, %f\n" % (centrox, centroy))
        f.write("%f\n\n" % radio)
            
        for punto in puntos:
            f.write("%f,%f\n" % (punto.x, punto.y))
        
        # Si no es el último cluster añadimos saltos de línea para diferenciarlos
        if(i < len(clusters) - 1):
            f.write("\n\n\n")
        
    f.close()  
    
    
    
# Extra: crear puntos a partir de circunferencias:
# Parámetros: nombre del archivo a crear, array con las circunferencias, número de puntos que queremos de cada una y rangos en los que queramos que puedan variar los puntos
# REF: https://gis.stackexchange.com/questions/76745/creating-a-circle-with-points
def generar_puntos_de_circunferencias(nombre_archivo, circunferencias, n_puntos_circunferencia, rango_x, rango_y):
    
    #Calculamos los puntos
    puntos = []
    
    for i in range(len(circunferencias)):
        x = circunferencias[i].centro.x
        y = circunferencias[i].centro.y
        radio = circunferencias[i].radio
        arc = (2 * math.pi) / n_puntos_circunferencia[i]
        
        for p in range(n_puntos_circunferencia[i]):
            px = (0*math.cos(arc * p)) - (radio*math.sin(arc * p))
            py = (radio*math.cos(arc * p)) + (0*math.sin(arc * p))
            px += x
            py += y
            puntos.append((px,py))
            
            
    # Los escribimos en un fichero      
    f = open("%s.csv" % nombre_archivo,"w+")
    
    if(len(puntos)/2>3):
        intervalo = rd.randint(3,math.floor(len(puntos)/2))
    else:
        intervalo = 2
        
    for i in range(len(puntos)):
        
        if(i%intervalo==0):
            x = puntos[i][0]
            y = puntos[i][1]
            x = x + rd.uniform(-rango_x, rango_x)
            y = y + rd.uniform(-rango_y, rango_y)
            
            f.write("%f,%f\n" % (x,y))
        else:
            f.write("%f,%f\n" % puntos[i])
    
    f.close() 
    

    
#########################################################################################################
############################################### Algoritmo ###############################################
#########################################################################################################



# Genera clusters iniciales aleatoriamente, y por cada uno de ellos realiza clustering. Se queda con los que tengan mayor grado de pertenencia. Por último elimina el posible ruido que haya podido quedar.
# Los parámetros están explicados en después del código
def clustering_circunferencias_con_incertidumbre(puntos, n_circunferencias, min_grado_pertenencia_clusters_finales, max_inicializaciones, razon_parada_clustering, max_iteraciones_clustering, min_grado_pertenencia_eliminar_punto_clustering, max_distancia_eliminar_punto_ruido, info):
    
    iterar = True
    n_iteraciones = 0
    
    # Guardamos el mejor resultado para que en caso de que se supere el máximo nº de iteraciones se devuelva este y no el último
    mejores_clusters = None
    mejores_clusters_min_grado = 0
    
    while(iterar and n_iteraciones < max_inicializaciones):
        
        # Inicializamos los clusters con circunferencias al azar a partir de todos los puntos
        clusters_iniciales = obtener_clusters_iniciales(puntos, n_circunferencias, min_grado_pertenencia_eliminar_punto_clustering, info)
        
        # Sobre esos clusters iniciales hacemos clustering para que se ajusten a los puntos
        clusters = clustering_circunferencias(puntos, clusters_iniciales, razon_parada_clustering, max_iteraciones_clustering, min_grado_pertenencia_eliminar_punto_clustering, info)
        
        # Si una de las circunferencias tiene menos de 3 puntos(los necesarios para calcular la nueva circunferencia), se sale de este bucle y se vuelve a comenzar.
        if(clusters == None):
            n_iteraciones = n_iteraciones + 1
            continue
        
        #Calculamos el mínimo grado de pertenencia de los puntos con su circunferencia asociada tras el clustering
        min_grado_actual = comprobar_grado_pertenencia(clusters)
        if(min_grado_actual == None):
            n_iteraciones = n_iteraciones + 1
            continue
        
        #Si se ha mejorado, se guarda el resultado (En la primera iteración siempre se guardará)
        if(mejores_clusters_min_grado < min_grado_actual):
            
            mejores_clusters = clusters
            mejores_clusters_min_grado = min_grado_actual
            
            # Si se ha alcanzado el menor grado de pertenencia de todos los puntos esperado por el usuario, se devuelve el resultado
            if(min_grado_actual > min_grado_pertenencia_clusters_finales):
                iterar = False
            
        n_iteraciones = n_iteraciones + 1
        
        
    #Info   
    if('n iteraciones inicializaciones' in info):
        cuadroLog.insert(INSERT, "Nº de iteraciones inicializaciones: " + str(n_iteraciones) + "\n")
        
    #Una vez que se hayan obtenido los mejores clusters posibles, se elimina el ruido, es decir, aquellos puntos cuya distancia a la circunferencia sea mayor que min_ruido
    if(max_distancia_eliminar_punto_ruido > 0):
        mejores_clusters = eliminar_ruido(mejores_clusters, max_distancia_eliminar_punto_ruido)
    
    # Info
    if('puntos eliminados por ruido' in info):
        puntos_clusters = []
        for cluster in mejores_clusters:
            puntos_clusters.extend(cluster.puntos)
            
        # Mostramos los puntos de la lista completa de puntos que no se encuentran en ningun cluster. Serán aquellos eliminados tanto en el clustering como en eliminar_ruido
        puntos_eliminados = []
        for punto in puntos:
            if(punto not in puntos_clusters):
                puntos_eliminados.append(punto)
                
        cuadroLog.insert(INSERT, "Puntos eliminados considerados como ruido: " + str(puntos_eliminados) + "\n")
        
    
    
    return mejores_clusters


# Funciones utilizadas por el algoritmo en orden de llamada:


# Crea n clusters iniciales al azar a partir de 3 puntos de la lista de todos los puntos.
# Parámetros: Array de todos los puntos, número de circunferencias totales, minimo grado de pertenencia clustering
def obtener_clusters_iniciales(puntos,n_circunferencias, min_grado_pertenencia_eliminar_punto_clustering, info): 
    
    # Primero se obtienen las circunferencias al azar
    circunferencias = []
    
    for i in range(n_circunferencias):
        
        puntos_aleatorios = rd.sample(puntos,3)
        circunferencia = obtener_circunferencia(puntos_aleatorios)
        circunferencias.append(circunferencia)
    
    # Se obtienen los clusters a partir de dichas circunferencias
    clusters_iniciales = obtener_clusters(circunferencias, puntos, min_grado_pertenencia_eliminar_punto_clustering, info)
    
    return clusters_iniciales


# Le asocia a cada circunferencia los puntos que tengan mayor grado de pertenencia con ella
# Parámetros: array de las circunferencias y todos los puntos
def obtener_clusters(circunferencias, puntos, min_grado_pertenencia_eliminar_punto_clustering, info):
    
    # Creamos clusters con las circunferencias pasadas
    clusters = []
    for circunferencia in circunferencias:
        clusters.append(Cluster(circunferencia))
    
    # Calculamos los grados de pertenencia de todos los puntos con todos los clusters
    for p in puntos:
        
        gradosPert = grados_pertenencia(p,circunferencias)
        
        # Si el máximo grado de pertenencia de un punto es menor que min_grado_pertenencia_clustering no se asigna a ningún cluster, y ya no se vuelve a considerar ese punto
        if(max(gradosPert) > min_grado_pertenencia_eliminar_punto_clustering):
            
            # Se le asigna al cluster con mayor grado de pertenencia
            indice = np.where(gradosPert == max(gradosPert))
            clusters[indice[0][0]].anyadir_punto(p)

    return clusters



# Devuelve un array con los grados de pertenencia de un punto a todas las circunferencias
# Parámetros: punto al que queremos calcular los grados de pertenencia, array de todas las circunferencias
def grados_pertenencia(p,circunferencias):
    
    #Calculamos la distancia del punto a cada circunferencia
    distancias = []
    
    for i in range(len(circunferencias)):
        d = np.sqrt(pow(p.x - circunferencias[i].centro.x,2) + pow(p.y - circunferencias[i].centro.y,2))
        distancias.append(abs(d - circunferencias[i].radio))
    
    #Calculamos el grado de pertencia a cada circunferencia con la siguiente fórmula: 100/distancia^2
    pertenencias = []
    
    for i in range(len(circunferencias)):
        if(distancias[i]==0.0):
            per = 100.0
        else:
            per = 100/(pow(distancias[i],2))
        pertenencias.append(per)
        
    #Normalizamos
    suma = sum(pertenencias)
    pertenencias = np.divide(pertenencias,suma)
    
    return pertenencias


# Obtiene una circunferencia mas cercana a partir de sus puntos asociados mediante clustering, pero estos pueden ir variando
# Parámetros: array con todos los puntos, los clusteres previos, la razon de parada, es decir, a partir de que diferencia minima entre circunferencias parará, y el máximo número de iteraciones que realizará.
def clustering_circunferencias(puntos, clusters_prev, razon_parada_clustering, max_iteraciones_clustering, min_grado_pertenencia_eliminar_punto_clustering, info):
    
    iterar = True
    iteraciones = 0;
    
    while(iterar and iteraciones <= max_iteraciones_clustering):
        
        diferencia = []
        
        # Obtenemos clusters mejores (Explicación en el método)
        clusters_post = obtener_siguientes_clusters(puntos, clusters_prev, min_grado_pertenencia_eliminar_punto_clustering, info)
        
        # Si una de las circunferencias tiene menos de 3 puntos(los necesarios para calcular la nueva circunferencia), no se realiza la clusterización, y se empezará con nuevas circunferencias aleatorias.
        if(clusters_post == None):
            return None;
        
        # Calculamos la diferencia entre los clusters previos y posteriores a partir de su centro y radio.
        for i in range(len(clusters_prev)):
            
            prev_centro = clusters_prev[i].circunferencia.centro
            prev_radio = clusters_prev[i].circunferencia.radio
            post_centro = clusters_post[i].circunferencia.centro
            post_radio = clusters_post[i].circunferencia.radio
            
            diferencia.append(np.sqrt((prev_centro.x - post_centro.x)**2 + (prev_centro.y - post_centro.y)**2) + abs(prev_radio - post_radio))
            
        # Si es menor que la razón, es el final de la clusterización
        if(max(diferencia) < razon_parada_clustering):    
            iterar = False
            
        # Si no ha parado, en la siguiente iteración el cluster previo será el calculado en esta
        clusters_prev = clusters_post
        iteraciones = iteraciones + 1
        
    #Info   
    if('n iteraciones clustering' in info):
        cuadroLog.insert(INSERT, "Nº de iteraciones clustering: " + str(iteraciones) + "\n")
        
        
    return clusters_post


# Se calculan los siguientes clusters a partir de los puntos asociados a cada uno de ellos.
def obtener_siguientes_clusters(puntos, clusters_prev, min_grado_pertenencia_eliminar_punto_clustering, info):
    
    nuevas_circunferencias = []
    
    for cluster_prev in clusters_prev:
        
        puntos_cluster = cluster_prev.puntos
        
        # Si una de las circunferencias tiene menos de 3 puntos(los necesarios para calcular la nueva circunferencia), no se realiza la clusterización.
        if(len(puntos_cluster) < 3):
            return None
        
        # Se obtiene la nueva circunferencia a partir de 3 puntos asociados a ella
        puntos_aleatorios = rd.sample(puntos_cluster,3)
        nuevas_circunferencias.append(obtener_circunferencia(puntos_aleatorios))
        
    # Se devuelve la nueva circunferencia con los nuevos puntos asociados a ella (Pueden ser los mismos que antes o no)
    clusters_post = obtener_clusters(nuevas_circunferencias, puntos, min_grado_pertenencia_eliminar_punto_clustering, info)
    
    return clusters_post



# Genera una circunferencia a partir de los 3 puntos pasados
# Parámetros: Array de 3 puntos
# REF: codewars.com/kata/give-the-center-and-the-radius-of-circumscribed-circle-a-warm-up-challenge/python
def obtener_circunferencia(puntos): 
    
    x1 = puntos[0].x
    x2 = puntos[1].x
    x3 = puntos[2].x
    y1 = puntos[0].y
    y2 = puntos[1].y
    y3 = puntos[2].y
    
    D = 2*(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
    
    # Comprobamos si los 3 puntos forman una recta, y si es así, quitamos los puntos iguales calculando un punto que los sustituya
    if(D==0):
        if(float(x1)==float(x2) and float(y1)==float(y2) and float(x1)==float(x3) and float(y1)==float(y3)):
            x2 = x1 + rd.random()
            y2 = y1 + rd.random()
            x3 = x2 - rd.random()
            y3 = y2 - rd.random()
            puntosNuevos = [puntos[0], Punto(x2,y2), Punto(x3,y3)]
        elif(float(x1)==float(x2) and float(y1)==float(y2)):
            distancia = np.sqrt(pow(x1-x3,2)+pow(y1-y3,2))
            x2 = x1 + distancia
            y2 = y1 + distancia
            puntosNuevos = [puntos[0], puntos[2], Punto(x2,y2)]
        elif(float(x1)==float(x3) and float(y1)==float(y3)):
            distancia = np.sqrt(pow(x1-x2,2)+pow(y1-y2,2))
            x3 = x1 + distancia
            y3 = y1 + distancia
            puntosNuevos = [puntos[0], puntos[1], Punto(x3,y3)]
        else:
            distancia = np.sqrt(pow(x2-x3,2)+pow(y2-y3,2))
            x1 = x2 + distancia 
            y1 = y2 + distancia
            puntosNuevos = [puntos[1], puntos[2], Punto(x1,y1)]
            
        
        return obtener_circunferencia(puntosNuevos)
    
    #Para hallar el centro (Ux, Uy), se calculan las mediatrices de los segmentos formados por los puntos y se escoge el valor donde corta  
    Ux = ((pow(x1,2) + pow(y1,2)) * (y2-y3) + (pow(x2,2) + pow(y2,2)) * (y3-y1) + (pow(x3,2) + pow(y3,2)) * (y1-y2))/D
    Uy = ((pow(x1,2) + pow(y1,2)) * (x3-x2) + (pow(x2,2) + pow(y2,2)) * (x1-x3) + (pow(x3,2) + pow(y3,2)) * (x2-x1))/D
    
    #Para hallar el diametro, calculamos las distancias Euclídeas entre puntos y el producto de ellas se divide entre D anteriormente definido
    AB = np.sqrt(float(pow(x2-x1,2) + pow(y2-y1,2)))
    BC = np.sqrt(float(pow(x3-x2,2) + pow(y3-y2,2)))
    AC = np.sqrt(float(pow(x3-x1,2) + pow(y3-y1,2)))
    diametro = (2*AB*BC*AC)/abs(D)
    
    circunferencia = Circunferencia(float(Ux),float(Uy),float(diametro/2))
    
    return circunferencia
    

# Devuelve el mínimo grado de pertenencia que tenga un punto de los puntos asociados a una circunferencia respecto a esa circunferencia de entre todas ellas.
def comprobar_grado_pertenencia(clusters):
    
    circunferencias = []
    grados_por_cluster = []
    grados_por_punto = []
    i = 0
    
    # Obtenemos todas las circunferencias
    for cluster in clusters:
        circunferencias.append(cluster.circunferencia)
        
    # Obtenemos el máximo grado de pertenencia de todos los puntos de un cluster a cada circunferencia
    for cluster in clusters:
        for p in cluster.puntos:
            grados_por_punto.append(max(grados_pertenencia(p,circunferencias)))
        
        # Si hay menos de 3 puntos, esta iteración se descarta y se vuelve a comenzar
        if(len(grados_por_punto) < 3):
            return None
        
        # De estos, obtenemos los mínimos grados  por cluster
        grados_por_cluster.append(min(grados_por_punto))
        
    # Por ultimo, devolvemos el menor de todos ellos
    return min(grados_por_cluster)



# Obtiene los nuevos clusters tras eliminar de ellos los puntos que se consideren como ruido
def eliminar_ruido(clusters, max_distancia_eliminar_punto_ruido):
    
    clusters_sin_ruido = []
    
    for cluster in clusters:
        clusters_sin_ruido.append(obtener_cluster_sin_ruido(cluster, max_distancia_eliminar_punto_ruido))
        
    return clusters_sin_ruido

# Dado un cluster, trata de obtener la circunferencia en las que la media de las distancias entre todos sus puntos y la circunferencia sea menor, y elimina los puntos que no cumplan cierto rango
def obtener_cluster_sin_ruido(cluster, max_distancia_eliminar_punto_ruido):
    
    
    clusters_finales = []
    puntos = cluster.puntos
    
    # Obtenemos todas las circunferencias posibles de recorrer los puntos 1 a 1
    # Por cada circunferencia:
    circunferencias = []
    for i in range(len(puntos) - 2):
    
        pts = [puntos[i],puntos[i+1],puntos[i+2]]
        circunferencias.append(obtener_circunferencia(pts))
    
    # Obtenemos todos los clusters posibles, pero elimiando de ellos aquellos puntos que sean considerados ruido
    for circunferencia in circunferencias:
        
        distancias = []
        centro = circunferencia.centro
        radio = circunferencia.radio
        
        # Calculamos las distancias de cada punto a la circunferencia
        for p in puntos:
            d = np.sqrt(pow(p.x - centro.x,2) + pow(p.y - centro.y,2))
            distancias.append(abs(d - radio))
            
        # Calculamos la media de esas distancias
        media = np.mean(distancias)
        
        # Seleccionamos los puntos cuya diferencia entre su distancia y la media sea menor que el min_ruido
        puntos_no_eliminados = []
        puntos_eliminados = []
        for i in range(len(distancias)):
            diferencia = abs(distancias[i] - media)
            if(diferencia < max_distancia_eliminar_punto_ruido):
                puntos_no_eliminados.append(puntos[i])
        
        # Si tuviesemos menos de 3 puntos no eliminados, volvemos a llamar al método con una mayor tolerancia para que coja mas puntos
        if(len(puntos_no_eliminados) < 3):
            return obtener_cluster_sin_ruido(cluster, max_distancia_eliminar_punto_ruido + 0.05)
        
        # Creamos un cluster con esa circunferenciua
        circun = Circunferencia(centro.x,centro.y,radio)
        clstr = Cluster(circun)
        clstr.anyadir_puntos(puntos_no_eliminados)
        clusters_finales.append(clstr)
        
    
    # Una vez que de un solo cluster tengamos muchos clusters sin ruidos, seleccionamos aquel cuya media de diferencias entre la distancia y la media sea menor
    max_media = float('inf')
    for i in range(len(clusters_finales)):
        
        centro = clusters_finales[i].circunferencia.centro
        radio = clusters_finales[i].circunferencia.radio
        puntos = clusters_finales[i].puntos
        distancias = []
        
        # Calculamos las distancias entre todos los puntos de un cluster y su circunferencia
        for p in puntos:
            d = np.sqrt(pow(p.x - centro.x,2) + pow(p.y - centro.y,2))
            distancias.append(abs(d - radio))
        
        # Calculamos la media
        media = np.mean(distancias)
        
        # Guardamos aquel cluster cuya media de distancias sea menor
        if(media < max_media):
            res = clusters_finales[i]
            max_media = media
    
    return res



##########################################################################################################################################
############################################### Interfaz Gráfica #########################################################################
##########################################################################################################################################


#Imports UI
from tkinter import *
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk
import tkinter.font as tkFont
from PIL import ImageTk,Image  
import os
import re

#Funciones
def verPuntos():
    if(localizacionPuntos.get()==""):
        messagebox.showinfo("Error...", "El campo de localización de archivo no puede estar vacío.")
    else:
        try:
            puntos = leer_puntos("%s" % localizacionPuntos.get())
        except:
            messagebox.showinfo("Error...", "El archivo no existe o no ha podido ser leido.")
            return 0
        try:
            dibujar_puntos(puntos)
            raiz3 = Tk()
            raiz3.title("Resultado")
            raiz3.iconbitmap("Data/Icono.ico")
            #Canvas
            canvas = Canvas(master = raiz3, width = 420, height = 280)  
            canvas.pack(fill="both", expand="True")
            #Imagen
            img = ImageTk.PhotoImage(Image.open("Data/puntosIniciales.png"), master=canvas)  
            canvas.create_image(0, 0, anchor=NW, image=img) 
            raiz3.mainloop()
        except:
            messagebox.showinfo("Error...", "El archivo no contiene los puntos correctamente, pruebe nuestro generador de puntos.")

def  ejecutar():
    cuadroLog.delete("1.0", END)
    cuadroLog.configure(state="normal")
    if(localizacionPuntos.get()=="" or cuadroNumCircunferencia.get()=="" or minGradoPertenencia.get()=="" or maxIteracionesAlgoritmo.get()=="" or razon.get()=="" or maxIteracionesClustering.get()=="" or minGradoPertenenciaClustering.get()=="" or minRuido.get()==""):
        messagebox.showinfo("Error...", "Alguno de los campos obligatorios está vacío.")
    else:
        try:
            puntos = leer_puntos("%s" % localizacionPuntos.get())
        except:
            messagebox.showinfo("Error...", "El archivo no existe o no ha podido ser leído.")
            return 0
        try:
            if(re.match("^[1-9][0-9]{0,5}$", cuadroNumCircunferencia.get()) is None or int(cuadroNumCircunferencia.get())<1):
                messagebox.showinfo("Error...", "El número de circunferencia tiene que ser un número mayor o igual a 1")
                return 0
            num_circunferencias = int(cuadroNumCircunferencia.get())
            if(re.match("^[0][.][0-9]{0,2}$", minGradoPertenencia.get()) is None or float(minGradoPertenencia.get())>=1.00 or float(minGradoPertenencia.get())<=0.00):
                messagebox.showinfo("Error...", "El mínimo grado de pertenencia tiene que ser mayor que 0.00 y menor que 1.00")
                return 0
            min_grado_pertenencia_clusters_finales = float(minGradoPertenencia.get())
            if(re.match("^[1-9][0-9]{0,6}$", maxIteracionesAlgoritmo.get()) is None or int(maxIteracionesAlgoritmo.get())<1):
                messagebox.showinfo("Error...", "El número máximo de inicializaciones tiene que ser un número mayor o igual a 1")
                return 0
            max_inicializaciones = int(maxIteracionesAlgoritmo.get())
            if(re.match("^[0-9]{1,3}[.][0-9]{0,2}$", razon.get()) is None or float(razon.get())>=1000.00 or float(razon.get())<=0.00):
                messagebox.showinfo("Error...", "La razón de parada tiene que ser mayor que 0.00 y menor que 1000.00")
                return 0
            razon_parada_clustering = float(razon.get())
            if(re.match("^[1-9][0-9]{0,6}$", maxIteracionesClustering.get()) is None or int(maxIteracionesClustering.get())<1):
                messagebox.showinfo("Error...", "El número máximo de iteraciones por clustering tiene que ser un número mayor o igual a 1")
                return 0
            max_iteraciones_clustering = int(maxIteracionesClustering.get())
            if(re.match("^[0][.][0-9]{0,2}$", minGradoPertenenciaClustering.get()) is None or float(minGradoPertenenciaClustering.get())>=1.00 or float(minGradoPertenenciaClustering.get())<=0.00):
                messagebox.showinfo("Error...", "El mínimo grado de pertenencia para eliminar un punto en el clustering tiene que ser mayor que 0.00 y menor que 1.00")
                return 0
            min_grado_pertenencia_eliminar_punto_clustering = float(minGradoPertenenciaClustering.get())
            if(re.match("^[\-]{0,1}[0-9]{1,5}[.]{0,1}[0-9]{0,2}$", minRuido.get()) is None or float(minRuido.get())>=10000.00 or float(minRuido.get())<=-10000.00):
                messagebox.showinfo("Error...", "La máxima distancia para eliminar un punto por ruido debe estar entre -10000.0 y 10000.0")
                return 0
            max_distancia_eliminar_punto_ruido = float(minRuido.get())
            info = []
            if(numIteIni.get()==1):
                info.append('n iteraciones inicializaciones')
            if(numIteClus.get()==1):
                info.append('n iteraciones clustering')
            if(punElim.get()==1):
                info.append('puntos eliminados por ruido')
            clusters = clustering_circunferencias_con_incertidumbre(puntos,num_circunferencias, min_grado_pertenencia_clusters_finales, max_inicializaciones, razon_parada_clustering,max_iteraciones_clustering, min_grado_pertenencia_eliminar_punto_clustering, max_distancia_eliminar_punto_ruido, info)
            cuadroLog.configure(state="disabled")
        except:
            messagebox.showinfo("Error...", "El algoritmo ha sufrido algún error, prueba de nuevo y si persiste contáctenos.")
            return 0
        if(localizacionPuntosResultado.get()!="" and nombrePuntosResultado.get()!=""):
            try:
                escribir_fichero_clusters(localizacionPuntosResultado.get()+"/"+nombrePuntosResultado.get(), clusters, 0)
                messagebox.showinfo("Archivo generado", "El archivo ha sido generado correctamente.")
            except:
                messagebox.showinfo("Error...", "La ruta especificada no es correcta.")
                return 0
        try:
            dibujar_clusters(clusters)
            raiz2 = Tk()
            raiz2.title("Resultado")
            raiz2.iconbitmap("Data/Icono.ico")
            #Canvas
            canvas = Canvas(master = raiz2, width = 420, height = 280)  
            canvas.pack(fill="both", expand="True")
            #Imagen
            img = ImageTk.PhotoImage(Image.open("Data/resultado.png"), master=canvas)  
            canvas.create_image(0, 0, anchor=NW, image=img) 
            raiz2.mainloop()
        except:
            messagebox.showinfo("Error...", "El clúster no ha podido ser dibujado, prueba de nuevo. En caso de persistir contáctenos.")
            return 0

def abrirArchivo():
    fichero=filedialog.askopenfilename(title="Abrir archivo de puntos", filetypes=(("Ficheros de CSV", "*.csv"),))
    localizacionPuntos.set(fichero)

def abrirCarpeta(localizacion):
    fichero=filedialog.askdirectory(title="Guardar en...")
    localizacion.set(fichero)
    
def abrirInstruccionesAnalizador():
    raiz2 = Tk()
    raiz2.title("Instrucciones del analizador de puntos")
    raiz2.iconbitmap("Data/Icono.ico")
    raiz2.geometry("550x650")
    #Fuentes
    fuenteTitulo = ("Helvetica",18,"bold")
    fuenteTexto = ("Helvetica",10)
    #Titulos
    instruccionTituloLabel = Label(raiz2, justify=LEFT, font=fuenteTitulo, text="Intrucciones:")
    instruccionTituloLabel.grid(row=0, column=0, sticky="w", padx=4, pady=4)
    #Texto
    instruccionTexto1Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="1º Seleccione la localización del archivo .CSV que desea abrir.")
    instruccionTexto1Label.grid(row=1, column=0, sticky="w", padx=4, pady=4)

    instruccionTexto2Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="2º Seleccione 'Ver puntos' si quiere ver los puntos del archivo escogido.")
    instruccionTexto2Label.grid(row=2, column=0, sticky="w", padx=4, pady=4)

    instruccionTexto3Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="3º Introduzca el número de circunferencias que contiene el archivo.")
    instruccionTexto3Label.grid(row=3, column=0, sticky="w", padx=4, pady=4)

    instruccionTexto4Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="4º Rellene el campo mínimo grado de pertenencia con un valor entre 0 y 1.")
    instruccionTexto4Label.grid(row=4, column=0, sticky="w", padx=4, pady=4)
    
    instruccionTexto5Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="5º Rellene el campo máximas inicializaciones con un valor mayor o igual que 1.")
    instruccionTexto5Label.grid(row=5, column=0, sticky="w", padx=4, pady=4)
    
    instruccionTexto6Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="6º Rellene el campo razón de parada de clastering con un valor mayor que 0.")
    instruccionTexto6Label.grid(row=6, column=0, sticky="w", padx=4, pady=4)
    
    instruccionTexto7Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="7º Rellene el campo máximo iteraciones de clustering con un valor mayor que 1")
    instruccionTexto7Label.grid(row=7, column=0, sticky="w", padx=4, pady=4)
    
    instruccionTexto8Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="8º Rellene el campo mínimo grado de pertenencia para eliminar punto de clustering\n    con un valor entre 0 y 1.")
    instruccionTexto8Label.grid(row=8, column=0, sticky="w", padx=4, pady=4)
    
    instruccionTexto9Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="9º Rellene el campo máxima distancia eliminar punto ruido, si es 0 o menor no se tendrá en cuenta.")
    instruccionTexto9Label.grid(row=9, column=0, sticky="w", padx=4, pady=4)

    instruccionTexto10Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="10º En caso de querer guardar el archivo de resultado, indique el nombre y\n     seleccione una ruta con el botón 'Abrir carpeta'.")
    instruccionTexto10Label.grid(row=10, column=0, sticky="w", padx=4, pady=4)
    
    instruccionTexto11Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="11º Marque las casillas en caso de querer imprimir parte de la información\n     que genera el algoritmo.")
    instruccionTexto11Label.grid(row=11, column=0, sticky="w", padx=4, pady=4)
    
    instruccionTexto12Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="12º En caso de tener dudas sobre los parámetros diríjase a la documentación\n     al apartado resultados.")
    instruccionTexto12Label.grid(row=12, column=0, sticky="w", padx=4, pady=4)
    
    instruccionTexto13Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="Parámetros recomendados:\n \n"+
                                    "num_circunferencias (K): Es muy importante que se ajuste a los datos.\n"+
                                    "min_grado_pertenencia_clusters_finales: casos sencillos: 0.8-0.9, casos complejos: 0.99.\n"+
                                    "max_inicializaciones: 10 - 200.\n"+
                                    "razon_parada_clustering: 0.01.\n"+
                                    "max_iteraciones_clustering: 3 - 30.\n"+
                                    "min_grado_pertenencia_eliminar_punto_clustering: 0.2.\n"+
                                    "max_distancia_eliminar_punto_ruido: 0.\n")
    instruccionTexto13Label.grid(row=13, column=0, sticky="w", padx=4, pady=4)
    
    
    raiz2.mainloop()
    
def abrirInstruccionesGenerador():
    raiz2 = Tk()
    raiz2.title("Instrucciones del analizador de puntos")
    raiz2.iconbitmap("Data/Icono.ico")
    raiz2.geometry("550x450")
    #Fuentes
    fuenteTitulo = ("Helvetica",18,"bold")
    fuenteTexto = ("Helvetica",10)
    #Titulos
    instruccionTituloLabel = Label(raiz2, justify=LEFT, font=fuenteTitulo, text="Intrucciones:")
    instruccionTituloLabel.grid(row=0, column=0, sticky="w", padx=4, pady=4)
    #Texto
    instruccionTexto1Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="1º Rellene el campo 'Circunferencia' siguiendo el patrón X.X, Y.Y, R.R donde X e Y son\n    las coordenadas del centro y R el radio de la circunferencia en decimal.")
    instruccionTexto1Label.grid(row=1, column=0, sticky="w", padx=4, pady=4)

    instruccionTexto2Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="2º Seleccione 'Añadir circunferencia' para añadirla a la lista inferior, en caso de\n    pulse el botón 'Limpiar' y vaciará la lista.")
    instruccionTexto2Label.grid(row=2, column=0, sticky="w", padx=4, pady=4)

    instruccionTexto3Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="3º Rellene el campo 'Número de puntos' para indicar cuandos puntos quiere de cada\n    circunferencia respectivamente.")
    instruccionTexto3Label.grid(row=3, column=0, sticky="w", padx=4, pady=4)

    instruccionTexto4Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="4º Seleccione 'Añadir número' para añadirla a la lista inferior, en caso de pulse el\n    botón 'Limpiar' y vaciará la lista.")
    instruccionTexto4Label.grid(row=4, column=0, sticky="w", padx=4, pady=4)

    instruccionTexto4Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="5º Rellene los campos 'Rango X' y 'Rango Y', para indicar el ruido máximo de los puntos\n    generados, en caso de no querer ruido dejar 0.0.")
    instruccionTexto4Label.grid(row=5, column=0, sticky="w", padx=4, pady=4)
    
    instruccionTexto5Label = Label(raiz2, justify=LEFT, font=fuenteTexto, text="6º Para guardar el archivo generado, indique el nombre y seleccione una ruta con el\n    botón 'Abrir carpeta'.")
    instruccionTexto5Label.grid(row=6, column=0, sticky="w", padx=4, pady=4)
    
    raiz2.mainloop()
    
def salir():
    respuesta = messagebox.askquestion("Salir", "¿Desea salir de la aplicación?")
    if(respuesta=="yes"):
        raiz.destroy()
        
def contacto():
    messagebox.showinfo("Contacto","En caso de duda o error en la aplicación contáctenos en:\n\nÁlvaro Aguilar Alhama: alvagualh@alum.us.es\nJosé Manuel Cobo Guerrero: joscobgue@alum.us.es")

#Raiz
raiz = Tk()
raiz.title("Clustering con incertidumbre")
raiz.geometry("700x700")
raiz.iconbitmap("Data/Icono.ico")

#Barra de menu
barraMenu=Menu(raiz)
raiz.config(menu=barraMenu)
archivoMenu=Menu(barraMenu, tearoff=0)
archivoMenu.add_command(label="Abrir archivo", command=abrirArchivo)
archivoMenu.add_separator()
archivoMenu.add_command(label="Instrucciones analizador", command=abrirInstruccionesAnalizador)
archivoMenu.add_command(label="Instrucciones generador", command=abrirInstruccionesGenerador)
archivoMenu.add_separator()
archivoMenu.add_command(label="Salir", command=salir)
barraMenu.add_cascade(label="Archivo", menu=archivoMenu)

ayudaMenu=Menu(barraMenu, tearoff=0)
ayudaMenu.add_command(label="Contacto", command=contacto)
barraMenu.add_cascade(label="Ayuda", menu=ayudaMenu)

#Notebook
notebook = ttk.Notebook(raiz)
notebook.pack()

#MiFrame
miFrame=Frame(notebook)
miFrame.pack(fill="both", expand="True")
notebook.add(miFrame, text="Obtención de circunferencias")

miFrame2=Frame(notebook)
miFrame2.pack(fill="both", expand="True")
notebook.add(miFrame2, text="Generador de puntos")

#######################################################Generador de circunferencias miFrame##############################################################

#Variables
localizacionPuntos = StringVar()
minGradoPertenencia = StringVar(raiz, value='0.95')
maxIteracionesAlgoritmo = StringVar(raiz, value='50')
razon = StringVar(raiz, value='0.01')
maxIteracionesClustering = StringVar(raiz, value='10')
minGradoPertenenciaClustering = StringVar(raiz, value='0.2')
minRuido = StringVar(raiz, value='0')
nombrePuntosResultado = StringVar()
localizacionPuntosResultado = StringVar()
numIteIni = IntVar()
numIteClus = IntVar()
punElim = IntVar()

#Cuadro de textos
cuadroLocalizacionPuntos = Entry(miFrame, textvariable=localizacionPuntos)
cuadroLocalizacionPuntos.grid(row=0, column=1)
localizacionPuntosLabel = Label(miFrame, text="Localizacion del archivo de puntos:")
localizacionPuntosLabel.grid(row=0, column=0, sticky="e", padx=4, pady=4)

cuadroNumCircunferencia = Spinbox(miFrame, from_=1, to=15)
cuadroNumCircunferencia.grid(row=2, column=1)
numCircunferenciaLabel = Label(miFrame, text="Número de circunferencias:")
numCircunferenciaLabel.grid(row=2, column=0, sticky="e", padx=4, pady=4)

cuadroMinGradoPertenencia = Spinbox(miFrame, from_=0, to=1, increment=0.01, format="%0.2f", textvariable=minGradoPertenencia)
cuadroMinGradoPertenencia.grid(row=3, column=1)
minGradoPertenenciaLabel = Label(miFrame, text="Mínimo grado de pertenencia de los clusters finales:")
minGradoPertenenciaLabel.grid(row=3, column=0, sticky="e", padx=4, pady=4)

cuadroMaxIteracionesAlgoritmo = Spinbox(miFrame, from_=1, to=2000, textvariable=maxIteracionesAlgoritmo)
cuadroMaxIteracionesAlgoritmo.grid(row=4, column=1)
maxIteracionesAlgoritmoLabel = Label(miFrame, text="Número máximo de inicializaciones aleatorias:")
maxIteracionesAlgoritmoLabel.grid(row=4, column=0, sticky="e", padx=4, pady=4)

cuadroRazon = Spinbox(miFrame, from_=0, to=100, increment=0.01, format="%0.2f", textvariable=razon)
cuadroRazon.grid(row=5, column=1)
razonLabel = Label(miFrame, text="Número de razón de parada de clustering:")
razonLabel.grid(row=5, column=0, sticky="e", padx=4, pady=4)

cuadroMaxIteracionesClustering = Spinbox(miFrame, from_=1, to=500, textvariable=maxIteracionesClustering)
cuadroMaxIteracionesClustering.grid(row=6, column=1)
maxIteracionesClusteringLabel = Label(miFrame, text="Número máximo de iteraciones del clustering:")
maxIteracionesClusteringLabel.grid(row=6, column=0, sticky="e", padx=4, pady=4)

cuadroMinGradoPertenenciaClustering = Spinbox(miFrame, from_=0, to=1, increment=0.01, format="%0.2f", textvariable=minGradoPertenenciaClustering)
cuadroMinGradoPertenenciaClustering.grid(row=7, column=1)
minGradoPertenenciaClusteringLabel = Label(miFrame, text="Número mínimo grado de pertenencia para eliminar un punto en el clustering:")
minGradoPertenenciaClusteringLabel.grid(row=7, column=0, sticky="e", padx=4, pady=4)

cuadroMinRuido = Spinbox(miFrame, from_=0, to=100000, increment=0.01, format="%0.2f", textvariable=minRuido)
cuadroMinRuido.grid(row=8, column=1)
minRuidoLabel = Label(miFrame, text="Número máxima distancia para eliminar un punto considerado como ruido:")
minRuidoLabel.grid(row=8, column=0, sticky="e", padx=4, pady=4)

cuadroNombrePuntosResultado = Entry(miFrame, textvariable=nombrePuntosResultado)
cuadroNombrePuntosResultado.grid(row=9, column=1)
nombrePuntosResultadoLabel = Label(miFrame, text="Nombre del archivo de resultado:")
nombrePuntosResultadoLabel.grid(row=9, column=0, sticky="e", padx=4, pady=4)

cuadroLocalizacionPuntosResultado = Entry(miFrame, textvariable=localizacionPuntosResultado)
cuadroLocalizacionPuntosResultado.grid(row=10, column=1)
localizacionPuntosResultadoLabel = Label(miFrame, text="Localizacion del archivo de resultado de puntos:")
localizacionPuntosResultadoLabel.grid(row=10, column=0, sticky="e", padx=4, pady=4)

boton1 = Checkbutton(miFrame, text="Mostrar número de iteraciones de cada inicialización", variable=numIteIni)
boton1.grid(row=13, column=0, sticky="w")
boton2 = Checkbutton(miFrame, text="Mostrar número de iteraciones por clustering", variable=numIteClus)
boton2.grid(row=14, column=0, sticky="w")
boton3 = Checkbutton(miFrame, text="Mostrar puntos eliminados por ruido", variable=punElim)
boton3.grid(row=15, column=0, sticky="w")

cuadroLog = Text(miFrame, width=80, height=8)
cuadroLog.grid(row=16, column=0, columnspan=3, sticky="e")
scrollVert=Scrollbar(miFrame, command=cuadroLog.yview)
scrollVert.grid(row=16, column=3, sticky="nsw")
cuadroLog.config(yscrollcommand=scrollVert.set)
cuadroLog.configure(state="disabled")

#Botones
botonVerPuntos = Button(miFrame, text="Ver puntos", command=lambda:verPuntos())
botonVerPuntos.grid(row=1,column=2,sticky="e",padx=2, pady=2)

botonAbrirArchivo = Button(miFrame, text="Abrir archivo", command=lambda:abrirArchivo())
botonAbrirArchivo.grid(row=1,column=1,sticky="w",padx=2,pady=2)

botonAbrirCarpeta = Button(miFrame, text="Abrir carpeta", command=lambda:abrirCarpeta(localizacionPuntosResultado))
botonAbrirCarpeta.grid(row=11,column=1,sticky="w",padx=2,pady=2)

botonMandar = Button(miFrame, text="Analizar", command=lambda:ejecutar())
botonMandar.grid(row=12,column=1,sticky="e",pady=4)

#########################################################Generador de puntos miFrame2####################################################################

#Funciones
def añadirCircunferencia():
    cuadroCircunferenciasAñadidas.configure(state="normal")
    if(circunferenciaPuntos.get()==""):
        messagebox.showinfo("Error...", "El campo circunferencia está vacío")
    else:
        regex = re.compile("^[\-]{0,1}[0-9]{1,4}[.]{0,1}[0-9]{0,2}[,]\s[\-]{0,1}[0-9]{1,4}[.]{0,1}[0-9]{0,2}[,]\s[0-9]{1,4}[.]{0,1}[0-9]{0,2}$")
        if re.match(regex, circunferenciaPuntos.get()) is None:
            messagebox.showinfo("Error...", "El campo circunferencia no cumple el patrón X.X, Y.Y, R.R")
        else:
            valores = circunferenciaPuntos.get().split(", ")
            imprimir = "(%.1f,%.1f),%.1f" % (float(valores[0]), float(valores[1]), float(valores[2]))
            cuadroCircunferenciasAñadidas.insert(INSERT, imprimir + "\n")
            cuadroCircunferenciasAñadidas.configure(state="disabled")
            circunferenciaPuntos.set("")
    
def añadirNumeroPuntos():
    cuadroNumeroPuntosAñadidos.configure(state="normal")
    if(circunferenciaNumPuntos.get()==""):
        messagebox.showinfo("Error...", "El campo número de puntos está vacío.")
    else:
        regex = re.compile("^[1-9][0-9]{0,4}$")
        if (re.match(regex,circunferenciaNumPuntos.get()) is None or int(circunferenciaNumPuntos.get())<3):
            messagebox.showinfo("Error...", "El campo número de puntos no es un número o es menor de 3.")
        else:
            imprimir = "%d" % int(circunferenciaNumPuntos.get())
            cuadroNumeroPuntosAñadidos.insert(INSERT, imprimir + "\n")
            cuadroNumeroPuntosAñadidos.configure(state="disabled")
            circunferenciaNumPuntos.set("")

def resetText(cuadro):
    cuadro.delete("1.0", END)

def generaPuntos():
    if(localizacionPuntosGenerados.get()=="" or nombrePuntosGenerados.get()==""):
        messagebox.showinfo("Error...", "La ruta o el nombre del archivo están vacíos.")
    else:
        nombreArchivo = localizacionPuntosGenerados.get() + "/" + nombrePuntosGenerados.get()
        listaCircunferencias = cuadroCircunferenciasAñadidas.get("1.0", END).split("\n")
        listaNumPuntos = cuadroNumeroPuntosAñadidos.get("1.0", END).split("\n")
        if(len(listaCircunferencias)!=len(listaNumPuntos)):
            messagebox.showinfo("Error...", "El número de circunferencias y número de puntos no coincide.")
        elif(len(listaCircunferencias)==2):
            messagebox.showinfo("Error...", "No has añadido circunferencias.")
        elif(len(listaNumPuntos)==2):
            messagebox.showinfo("Error...", "No has añadido número de puntos a generar.")
        else:
            circunferencias = []
            nPuntos = []
            for i in range(len(listaCircunferencias) - 2):
                regexCircun = re.compile("^[(][\-]{0,1}[0-9]{1,4}[.]{0,1}[0-9]{0,2}[,][\-]{0,1}[0-9]{1,4}[.]{0,1}[0-9]{0,2}[)][,][0-9]{1,4}[.]{0,1}[0-9]{0,2}$")
                if re.match(regexCircun, listaCircunferencias[i]) is None:
                    messagebox.showinfo("Error...", "Hay una circunferencia en puntos añadidos que no cumple el formato dado.")
                    return 0
                trozos1 = listaCircunferencias[i].split("),")
                trozos2 = trozos1[0].replace("(","").split(",")
                circunferencias.append(Circunferencia(float(trozos2[0]), float(trozos2[1]), float(trozos1[1])))
                regexNum = re.compile("^[1-9][0-9]{0,4}$")
                if (re.match(regexNum, listaNumPuntos[i]) is None or int(listaNumPuntos[i])<3):
                    messagebox.showinfo("Error...", "Hay un número de puntos añadidos que no es un número o es menor que 3.")
                    return 0
                nPuntos.append(int(listaNumPuntos[i]))
            try:
                generar_puntos_de_circunferencias(nombreArchivo, circunferencias, nPuntos, float(rangoX.get()), float(rangoY.get()))
                localizacionPuntos.set(nombreArchivo + ".csv")
                messagebox.showinfo("Generado correctamente", "El fichero de puntos " + nombrePuntosGenerados.get() + " se ha generado correctamente.")
            except:
                messagebox.showinfo("Error...", "Se ha producido un error al generar los puntos, pruebe de nuevo. En caso de persistir contáctenos.")
                
#Variables
circunferenciaPuntos = StringVar(raiz, value="0.0, 0.0, 1.0")
nombrePuntosGenerados = StringVar()
localizacionPuntosGenerados = StringVar()
circunferenciaNumPuntos = StringVar()
rangoX = StringVar()
rangoY = StringVar()
    
#Cuadro de textos
cuadroCircunferenciaPuntos = Entry(miFrame2, textvariable=circunferenciaPuntos)
cuadroCircunferenciaPuntos.grid(row=0, column=1)
circunferenciaPuntosLabel = Label(miFrame2, text="Circunferencia:")
circunferenciaPuntosLabel.grid(row=0, column=0, sticky="e", padx=4, pady=4)

cuadroCircunferenciasAñadidas = Text(miFrame2, width=16, height=5)
cuadroCircunferenciasAñadidas.grid(row=2, column=1, sticky="e")
circunferenciasAñadidasLabel = Label(miFrame2, text="Circunferencias añadidas:")
circunferenciasAñadidasLabel.grid(row=2, column=0, sticky="e", padx=4, pady=4)
scrollVert=Scrollbar(miFrame2, command=cuadroCircunferenciasAñadidas.yview)
scrollVert.grid(row=2, column=2, sticky="nsw")
cuadroCircunferenciasAñadidas.config(yscrollcommand=scrollVert.set)
cuadroCircunferenciasAñadidas.configure(state="disabled")

cuadroCircunferenciaNumPuntos = Spinbox(miFrame2, from_=3, to=1000, textvariable=circunferenciaNumPuntos)
cuadroCircunferenciaNumPuntos.grid(row=4, column=1)
circunferenciaNumPuntosLabel = Label(miFrame2, text="Número de puntos:")
circunferenciaNumPuntosLabel.grid(row=4, column=0, sticky="e", padx=4, pady=4)

cuadroNumeroPuntosAñadidos = Text(miFrame2, width=16, height=5)
cuadroNumeroPuntosAñadidos.grid(row=6, column=1, sticky="e")
numeroPuntosAñadidosLabel = Label(miFrame2, text="Número de puntos añadidos:")
numeroPuntosAñadidosLabel.grid(row=6, column=0, sticky="e", padx=4, pady=4)
scrollVert=Scrollbar(miFrame2, command=cuadroNumeroPuntosAñadidos.yview)
scrollVert.grid(row=6, column=2, sticky="nsw")
cuadroNumeroPuntosAñadidos.config(yscrollcommand=scrollVert.set)
cuadroNumeroPuntosAñadidos.configure(state="disabled")

cuadroRangoX = Spinbox(miFrame2, from_=0, to=100, increment=0.1, format="%.1f", textvariable=rangoX)
cuadroRangoX.grid(row=8, column=1)
rangoXLabel = Label(miFrame2, text="Rango de variación de X:")
rangoXLabel.grid(row=8, column=0, sticky="e", padx=4, pady=4)

cuadroRangoY = Spinbox(miFrame2, from_=0, to=100, increment=0.1, format="%.1f", textvariable=rangoY)
cuadroRangoY.grid(row=9, column=1)
rangoYLabel = Label(miFrame2, text="Rango de variación de Y:")
rangoYLabel.grid(row=9, column=0, sticky="e", padx=4, pady=4)

cuadroNombrePuntosGenerados = Entry(miFrame2, textvariable=nombrePuntosGenerados)
cuadroNombrePuntosGenerados.grid(row=10, column=1)
nombrePuntosGeneradosLabel = Label(miFrame2, text="Nombre del archivo de puntos generados:")
nombrePuntosGeneradosLabel.grid(row=10, column=0, sticky="e", padx=4, pady=4)

cuadroLocalizacionPuntosGenerados = Entry(miFrame2, textvariable=localizacionPuntosGenerados)
cuadroLocalizacionPuntosGenerados.grid(row=11, column=1)
localizacionPuntosGeneradosLabel = Label(miFrame2, text="Localizacion del archivo de puntos generados:")
localizacionPuntosGeneradosLabel.grid(row=11, column=0, sticky="e", padx=4, pady=4)

#Botones
botonAñadirCircunferencia = Button(miFrame2, text="Añadir circunferencia", command=lambda:añadirCircunferencia())
botonAñadirCircunferencia.grid(row=1,column=2,sticky="w",padx=2, pady=2)

botonLimpiarCircunferencia = Button(miFrame2, text="Limpiar", command=lambda:resetText(cuadroCircunferenciasAñadidas))
botonLimpiarCircunferencia.grid(row=3,column=2,sticky="w",padx=2, pady=2)

botonAñadirNumPuntos = Button(miFrame2, text="Añadir número", command=lambda:añadirNumeroPuntos())
botonAñadirNumPuntos.grid(row=5,column=2,sticky="w",padx=2, pady=2)

botonLimpiarNumPuntos = Button(miFrame2, text="Limpiar", command=lambda:resetText(cuadroNumeroPuntosAñadidos))
botonLimpiarNumPuntos.grid(row=7,column=2,sticky="w",padx=2, pady=2)

botonAbrirCarpeta = Button(miFrame2, text="Abrir carpeta", command=lambda:abrirCarpeta(localizacionPuntosGenerados))
botonAbrirCarpeta.grid(row=12,column=1,sticky="w",padx=2,pady=2)

botonGeneraPuntos = Button(miFrame2, text="Generar puntos", command=lambda:generaPuntos())
botonGeneraPuntos.grid(row=13,column=1,sticky="w",padx=2,pady=6)

raiz.mainloop()


# ## Funcionamiento:
# 
# 
# > **clustering_circunferencias_con_incertidumbre**(puntos, n_circunferencias, min_grado_pertenencia_clusters_finales, max_inicializaciones, razon_parada_clustering, max_iteraciones_clustering, min_grado_pertenencia_eliminar_punto_clustering,    max_distancia_eliminar_punto_ruido, info)
# 
# &nbsp;
# 
#  - **puntos**: array de los puntos sobre los que se quiere hacer clustering.
# 
#  - **num_circunferencias (K)**: nº estimado por el usuario de la cantidad de clusters que hay.
# 
#  - **min_grado_pertenencia_clusters_finales**: No se harán más inicializaciones aleatorias si el mínimo grado de pertenencia de todos los puntos de los clusters al que estan asociados es mayor que este valor. Rango: (0-1).
# 
#  - **max_inicializaciones**: Máximo nº de veces que se inicializará aleatoriamente las circunferencias iniciales.
# 
#  - **razon_parada_clustering**: Si la diferencia entre un cluster y el siguiente que se calcule es menor que este valor, para la clusterización. Rango: (0 - inf).
# 
#  - **max_iteraciones_clustering**: Nº máximo de iteraciones de clustering por cada inicialización.
# 
#  - **min_grado_pertenencia_eliminar_punto_clustering**: Si el mínimo grado de pertenencia de un punto con todas las circunferencias es menor que este valor, no se añade a ningun cluster. Rango: (0-1) Cuanto mayor sea, más probable es que se eliminen puntos que se encuentren entre varias circunferencias.
# 
#  - **max_distancia_eliminar_punto_ruido**: distancia máxima de un punto a su cluster para que no sea considerado ruido y se elimine. Rango:(0, +inf) Cuanto menor sea, más probable es que se elimine los puntos más alejados a un cluster. Si se le pasa un valor de 0 o menor, no se ejecutará la función eliminar_ruido().
#  
#  &nbsp;
#  
# **Parámetros recomendados:**
#  
#  - **num_circunferencias (K)**: Es muy importante que se ajuste a los datos, ya que, de lo contrario, los resultados serán erróneos.
#  
#  - **min_grado_pertenencia_clusters_finales**: En casos sencillos es válido un valor en torno a 0.8-0.9 para aumentar la velocidad, en cambio, en casos complejos se recomienda subir hasta 0.99.
#  
#  - **max_inicializaciones**: Al igual que con el parámetro anterior, en casos sencillos, con unas 10 será suficiente. A mayor complejidad, cuanto más alto sea el valor, mayor seguridad habrá de encontrar un buen resultado, encontrando casos en lo que lo recomendable sean 200, teniendo en cuenta que esto tendrá un coste superior en tiempo.
#  
#  - **razon_parada_clustering**: Cuanto más pequeño sea, se tendrá más seguridad de que se ha llegado a una solución óptima y que no variaran más las circunferencias, pero más tardará en parar. Es por ello que con un valor de 0.01 se considera adecuado.
#  
#  - **max_iteraciones_clustering**: Dependiendo de la complejidad, se recomienda un valor de entre 3 y 30.  A mayor valor más tardará en ejecutarse.
#  
#  - **min_grado_pertenencia_eliminar_punto_clustering**: Este es un caso complejo, y que se debe ajustar a los puntos iniciales. Si es un caso sencillo, es mejor dejar un valor pequeño, de 0.2 por ejemplo, ya que no influirá mucho. Pero si se detecta ruido entre las circunferencias, se puede aumentar hasta que este desaparezca. Debe tenerse en cuenta que dependerá de K, ya que, a más circunferencias, probablemente menor sea los grados de pertenencia de los puntos
# 
# - **max_distancia_eliminar_punto_ruido**: Se recomienda dejarlo en 0 o menos para que no se tenga en cuenta, pero si se detecta que hay ruido, se debe buscar un valor que se aproxime a la distancia de ese ruido al clúster al que pertenece. Utilizar esta función tiene un impacto importante en el rendimiento.

# In[ ]:




