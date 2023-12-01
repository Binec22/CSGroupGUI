from math import cos,sin,radians
import numpy as np
from pyproj import Transformer
from pyproj import Proj
import folium
import os
import random
import xml.etree.ElementTree as et
from folium import plugins
from datetime import datetime
from folium.plugins import TimestampedGeoJson


#####################################################
### Fonctions qui permettent de faire les calculs ###
#####################################################

def getAltitude(z0,vz0,list_temps,type_alt = 'linéaire',orientation = 'ascendant'):
    """
    ENTREES :
        z0 = altitude initiale
        vz0 = norme de la vitesse en altitude
        list_temps = liste des temsp calculée par les fonctions de trajectoires
        type = type de changement en altitude désiré
            - 'linéaire' = montée ou descente linéaire
            - 'parabole' = montée ou descente suivant l'equation d'une trajectoire parabolique (définir les paramètres d'accélération)
        orientation = 'ascendant' ou 'descendant'

    SORTIES :
        liste d'altitude en z de la trajectoire
    """

    # determination du coefficient de la vitesse du changement d'altitude
    vz0 = vz0 - 2*(orientation != 'ascendant')*vz0

    if type_alt == 'linéaire':
        return [(z0 + vz0*dt)*(z0 + vz0*dt >=0) + 0 for dt in list_temps]
    
    elif type_alt == 'parabole':
        # prise en compte ou non de l'accélération de la pesanteur
        str_pes = str(input("Trajectoire soumise seulement à la pesanteur ? (Oui/Non): "))

        if str_pes == 'Oui':
            return [(z0 - 0.5*9.81*dt**2 + vz0*dt)*(z0 - 0.5*9.81*dt**2 + vz0*dt >= 0) + 0 for dt in list_temps]
        
        elif str_pes == 'Non':
            # determination du paramètre de la parabole par l'utilisateur
            p = float(input("Paramètre de la parabole (p >= 0): "))

            if p < 0:
                raise(ValueError("Le paramètre de la parabole doit être positif ou nul"))
            
            # détermination de l'orientation de la trajectoire
            para_or = str(input("Orientation de la parabole ? (Haut/Bas): "))
            print(para_or)

            p = p - 2*(para_or == 'Bas')*p
            
            return [(z0 + p*dt**2 + vz0*dt)*(z0 + p*dt**2 + vz0*dt >= 0) + 0 for dt in list_temps] 



def rectiligne(x0,y0,v0,cap,temps_acquisition,pas_de_temps):
   """
   v0 = vitesse initiale, en m/s
   pas_de_temps = pas de temps entre deux points successifs, en s
   temps_aquisition = temps total de l'acquisition, en s

   SORTIE : liste des positions en x
            liste des positions en y
            liste des temps 
            liste des caps à tout instant
   """

   # transformation du cap en angle par rapport à l'axe x parallèle à l'axe Ouest-Est
   # conversion en radians
   phi = (-cap+90)*np.pi/180

   # initialisation du vecteur vitesse initiale
   v = v0*np.array([cos(phi),sin(phi)]).reshape(-1,1)

   # calcul de la matrice des temps
   dts = np.arange(0, temps_acquisition + pas_de_temps , pas_de_temps).reshape(1,-1)

   # calcul de la matrice des positions
   vecteurs = v@dts

   # mise en forme de la matrice des positions
   # ajout de la position à t = 0
   depart = np.array([x0,y0]).reshape(-1,1)
   ans = depart + vecteurs

   list_x , list_y = list(ans[0]) , list(ans[1])
   list_cap = [cap for _ in range(len(list_x))]

   return  list_x , list_y , dts.tolist()[0] , list_cap



def circulaire(x0,y0,v0,cap,rayon_courbure,temps_acquisition,pas_de_temps,orientation_virage):
   """
   v0 = vitesse initiale, en m/s
   rayon_courbure = rayon du cercle de la trajectorie, en m
   pas_de_temps = pas de temps entre deux points successifs, en s
   temps_aquisition = temps total de l'acquisition, en s

   SORTIE : liste des positions en x
            liste des positions en y
            liste des temps
            liste des caps à tout instant
   """
   R = rayon_courbure

   # transformation du cap en angle par rapport à l'axe x parallèle à l'axe Ouest-Est
   # convertion en radians
   phi=(-cap+90)*np.pi/180

   # definition du sens de parcours du cercle 
   if orientation_virage == 'Droite':
      R = -R
   else :
       R = R
   
   # calcul de la matrice des temps
   dts = np.arange(0 , temps_acquisition+pas_de_temps , pas_de_temps).reshape(1,-1)

   # initialisation du vecteur position
   traj=phi+v0*np.repeat(dts,2,axis=0)/R

   # calcul de la matrice de position
   traj=np.concatenate((x0+R*np.sin(traj[0,:]).reshape(1,-1)-R*sin(phi), y0-R*np.cos(traj[1,:]).reshape(1,-1)+R*cos(phi)),axis=0)

   dcap = v0/(R*pas_de_temps)
   
   list_cap = [(cap - i*dcap)%360 for i in range(len(traj[0]))]
   
   return list(traj[0]) , list(traj[1]) , dts.tolist()[0] , list_cap



def zigzag(x0,y0,v0,cap,temps_acquisition,pas_de_temps,nombre_lignes,longueur_ligne,angle_changement,virage_1,type_zigzag):
    """
    fonction zigzag unique, définition de paramètres suivant le type de zigzag voulu
    v0 = vitesse initiale, en m/s
    pas_de_temps = pas de temps entre deux points successifs, en s
    temps_aquisition = temps total de l'acquisition, en s
    nombre_lignes = le nombre de lignes que fera le zigzag si ce paramètre a été choisi
    longueur_ligne = longueur d'une ligne si ce paramètre a été choisi
    angle_changement = angle entre chaque zigzag

    type_zigzag = type de zigzag voulu
        - 'n_adef', pour un zigzag avec un angle de changement défini et un nombre de lignes défini
        - 'n_alea', pour un zigzag avec un angle de changement aléatoire mais un nombre de lignes défini
        - 'l_adef', pour un zigzag avec un angle de changement défini et une longueur de portions droites définies
        - 'l_alea', pour un zigzag avec un angle de changement aléatoire mais une longueur de portions droites définies

    SORTIE : liste des positions en x
             liste des positions en y
             liste des temps
             liste des caps à tout instant
    """

    Temps = [0]
    list_cap = [cap]

    # cas nombre de lignes défini, angle de changement défini
    if type_zigzag == 'n_adef':
        # calcul du temps de parcours d'une portion droite de la trajectoire
        temps_rect = temps_acquisition/nombre_lignes

        # creation d'une liste contenant les deux caps possibles du bateau suivant le zigzag #
        diff_cap = [cap,cap + 180 + (2*(virage_1 == 'Gauche') - 1)*angle_changement]

        if diff_cap[1] >= 360:
            diff_cap[1] -= 360
   
        # initialisation des listes des positions en x et y
        listePosition_x , listePosition_y = [x0] , [y0]

        # initialisation du cap utilisé la portion droite étudiée
        cap_used = cap


    # cas longueur de lignes définie, angle de changement défini
    if type_zigzag == 'l_adef':
        # calcul du temps de parcours d'une portion droite de longueur l_lignes
        temps_rect = longueur_ligne/v0

        if temps_rect > temps_acquisition:
            raise ValueError("Temps d'acquisition trop petit : augmenter temps d'acquisition ou diminuer temps parcours des lignes")
        
        else:
            # calcul du nombre de portions droites complètes de la trajectoire
            nombre_lignes = int(temps_acquisition//temps_rect)

            # creation d'une liste contenant les deux caps possibles du bateau suivant le zigzag #
            diff_cap = [cap,cap + 180 + (2*(virage_1 == 'Gauche') - 1)*angle_changement]

            if diff_cap[1] >= 360:
                diff_cap[1] -= 360
            
            # initialisation des listes des positions en x et y
            listePosition_x , listePosition_y = [x0] , [y0]

            # initialisation du cap utilisé la portion droite étudiée
            cap_used = cap


    # cas nombre de lignes défini, angle de changement aléatoire
    if type_zigzag == 'n_alea':
        # calcul du temps de parcours d'une portion droite de la trajectoire
        temps_rect = temps_acquisition/nombre_lignes      

        # initialisation des listes des positions en x et y
        listePosition_x , listePosition_y = [x0] , [y0]

        # initialisation du cap utilisé la portion droite étudiée
        cap_used = cap

        # initialisation du compteur d'orientation du virage
        cpt_virage = random.choice([-1,1])


    # cas longueur de lignes définie, angle de changement aléatoire
    if type_zigzag == 'l_alea':
        # calcul du temps de parcours d'une portion droite de longueur l_lignes
        temps_rect = longueur_ligne/v0

        # calcul du nombre de portions droite complète de la trajectoire
        nombre_lignes = int(temps_acquisition//temps_rect)      

        # initialisation des listes des positions en x et y
        listePosition_x , listePosition_y = [x0] , [y0]

        # initialisation du cap utilisé la portion droite étudiée
        cap_used = cap

        # initialisation du compteur d'orientation du virage
        cpt_virage = random.choice([-1,1])


    #----------------MAIN----------------#
    # boucle de calcul de la trajectoire #
    #------------------------------------#
    
    for _ in range(0,nombre_lignes):
        # parcours d'une portion droite avec un cap 'cap_used'
        x,y,dts,lc = rectiligne(listePosition_x[-1],listePosition_y[-1],v0,cap_used,temps_rect-pas_de_temps,pas_de_temps)
        x.pop(0)
        y.pop(0)
        dts.pop(0)
        lc.pop(0)

        # concatenantion des positions en x et y
        listePosition_x += x
        listePosition_y += y
        list_cap += lc

        # ajout des pas de temps de la ligne parcourue dans la liste des temps
        dts2 = [t + Temps[-1] + pas_de_temps for t in dts]
        Temps += dts2

        if type_zigzag == 'n_adef' or type_zigzag == 'l_adef':
            # changement de cap pour la ligne suivante
            if cap_used == diff_cap[0]:
                cap_used = diff_cap[1]
            else:
                cap_used = diff_cap[0]
        
        elif type_zigzag == 'n_alea' or type_zigzag == 'l_alea':
            # definition aleatoire de l'angle de changement de cap
            alp = random.uniform(0,360)

            # conversion de l'angle de changement en un nouveau cap
            cap_used = cap + 180 + cpt_virage*alp
            if cap_used >= 360:
                cap_used -= 360
            
            # mise à jour du compteur d'orientation des virages
            cpt_virage = -1*(cpt_virage == 1) + 1*(cpt_virage == -1)
    

    if type_zigzag == 'l_adef' or type_zigzag == 'l_alea':
        if temps_acquisition%temps_rect != 0:
            # calcul du temps réstant après parcours des lignes complètes 
            t_last_line = temps_acquisition - nombre_lignes*temps_rect

            # parcours de la dernière portion de ligne droite
            x,y,dts,lc = rectiligne(listePosition_x[-1],listePosition_y[-1],v0,cap_used,t_last_line,pas_de_temps)
            x.pop(0)
            y.pop(0)
            dts.pop(0)
            lc.pop(0)

            # concatenantion des positions en x et y de la dernière portion droite
            listePosition_x += x
            listePosition_y += y
            list_cap += lc

            # ajout des pas de temps de la dernière portion droite dans la liste des temps
            dts2 = [t + Temps[-1] + pas_de_temps for t in dts]
            Temps += dts2
    
    return listePosition_x,listePosition_y,Temps,list_cap



def zigzag_full_alea(x0,y0,v0,cap,temps_acquisition,pdt):
   """
   fonction zigzag à parcours 100% aléatoire

   v0 = vitesse initiale, en m/s
   pdt = pas de temps entre deux points successifs, en s
   temps_aquisition = temps total de l'acquisition, en s

   SORTIE : liste des positions en x
            liste des positions en y
            lsite des temps   
   """
   Temps = [0]
   temps_restant = temps_acquisition

   # initialisation du premier temps de parcours
   temps_rect = random.uniform(1,temps_restant/4)

   # temps restant après parcours de la première ligne droite
   temps_restant = temps_restant - temps_rect

   # initialisation des listes des positions en x et y
   listePosition_x , listePosition_y , list_cap = [x0] , [y0] ,[cap]

   # initialisation du cap utilisé la portion droite étudiée
   cap_used = cap

   # initialisation aleatoire du compteur d'orientation du virage
   cpt_virage = random.choice([-1,1])

   #initialisation du compteur de lignes effectuées
   cpt_lignes = 0

   while temps_restant > 0:
      # mise à jour du compteur de lignes
      cpt_lignes += 1

      # parcours d'une portion droite du zigzag
      x,y,dts,lc = rectiligne(listePosition_x[-1],listePosition_y[-1],v0,cap_used,temps_rect,pdt)
      x.pop(0)
      y.pop(0)
      dts.pop(0)
      lc.pop(0)
    
      # concatenation des positions en x et y
      listePosition_x += x
      listePosition_y += y
      list_cap += lc

      # ajout des pas de temps de la ligne parcourue dans la liste des temps
      dts2 = [t + Temps[-1] + pdt for t in dts]
      Temps += dts2

      # definition aleatoire de l'angle de changement de cap
      alp = random.uniform(0,360)

      # conversion de l'angle de changement en un nouveau cap
      cap_used = cap + 180 + cpt_virage*alp
      if cap_used >= 360:
         cap_used -= 360

      # mise à jour du compteur d'orientation des virages
      cpt_virage = -1*(cpt_virage == 1) + 1*(cpt_virage == -1)

      # mise à jour du temps de parcours de la prochaine ligne droite
      temps_rect = random.uniform(1,temps_restant)

      # mise à jour du temps restant après parcours de la prochaine ligne droite
      temps_restant = temps_restant - temps_rect
   
   return listePosition_x,listePosition_y,Temps,list_cap



def trajectoireArc(x0,y0,v0,cap,rayon_courbure,temps_acquisition,pas_de_temps,angle_rotation,orientation_virage):
   """ 
   v0 en m/s
   rayon_courbure en metres
   temps_acquisition en s
   pas_de_temps en s
   angle_rotation en degres

   SORTIE : liste des positions en x
            liste des positions en y
            liste des temps
            liste des caps à tout instant
   """
   
   R = rayon_courbure

   # calcul de la longueur l'arc de cercle
   d = radians(angle_rotation)*R

   # calcul du temps de parcours de l'arc de cercle
   temps_parcours = d/v0

   if temps_parcours > temps_acquisition:
      return circulaire(x0,y0,v0,cap,rayon_courbure,temps_acquisition,pas_de_temps,orientation_virage)

   else:
      # parcours de l'arc de cercle
      list_x,list_y,list_t,list_cap =  circulaire(x0,y0,v0,cap,rayon_courbure,temps_parcours,pas_de_temps,orientation_virage)

      # calcul du cap après parcours de l'arc de cercle
      if orientation_virage != 'Gauche':
         new_cap = cap + angle_rotation
      else:
         new_cap = cap - angle_rotation

      # calcul du temps restant après parcours de l'arc
      temps_restant = temps_acquisition - temps_parcours - pas_de_temps

      # parcours d'une ligne droite après l'arc de cercle
      x,y,dt,lc = rectiligne(list_x[-1],list_y[-1],v0,new_cap,temps_restant,pas_de_temps)

      x.pop(0)
      y.pop(0)
      dt.pop(0)
      lc.pop(0)

      list_x += x
      list_y += y
      list_cap += lc

      # ajout des pas de temps de la seconde ligne droite dans la liste de temps
      dt2 = [t + list_t[-1] for t in dt]
      list_t += dt2

      return list_x , list_y , list_t , list_cap



def calculArc(x0,y0,v0,cap,rayon_courbure,pas_de_temps,angle_rotation,orientation_virage = 'gauche'):
   """ 
   Fonction ne prenant pas en compte le temps d'acquisition, utilisée dans les fonctions Hippodrome et trajectoireU
   v0 en m/s
   rayon_courbure en metres
   pas_de_temps en s
   angle_rotation en degres

   SORTIE : liste des positions en x
            liste des positions en y
            liste des temps
            liste des caps à tout instant
   """
   
   R = rayon_courbure

   # calcul de la longueur l'arc de cercle
   d = radians(angle_rotation)*R

   # calcul du temps de parcours de l'arc de cercle
   temps_parcours = d/v0

   return circulaire(x0,y0,v0,cap,rayon_courbure,temps_parcours-pas_de_temps,pas_de_temps,orientation_virage)



def hippodrome(x0,y0,v0,cap,rayon_courbure,temps_acquisition,pas_de_temps,orientation_virage):
   """
   v0 = vitesse initiale, en m/s
   rayon_courbure = rayon du demi-cercle de la trajectorie, en m
   pas_de_temps
   pas_de_temps = pas de temps entre deux points successifs, en s
   temps_aquisition = temps total de l'acquisition, en s

   SORTIE : liste des positions en x
            liste des positions en y
            liste des temps
            liste des caps à tout instant
   """
   
   # calcul de la longueur du demi-cercle
   d = radians(180)*rayon_courbure

   # calcul du temps de parcours du demi cercle
   tps_arc = d/v0

   if tps_arc > temps_acquisition:
      raise(ValueError("Temps d'acquisition trop petit"))
   
   else:
      # calcul du temps restant 
      tps_restant = temps_acquisition - tps_arc

      # definition aléatoire du temps de parcours des parties rectiligne de la trajectoire
      tps_rect1 = random.uniform(pas_de_temps,tps_restant)
      tps_rect2 = tps_restant - tps_rect1

      # parcours de la première ligne droite
      list_x , list_y , Temps , list_cap = rectiligne(x0,y0,v0,cap,tps_rect1,pas_de_temps)

      # parcours du demi cercle
      la_x,la_y,dts,lc = calculArc(list_x[-1],list_y[-1],v0,cap,rayon_courbure,pas_de_temps,180,orientation_virage)
      la_x.pop(0)
      la_y.pop(0)
      dts.pop(0)
      lc.pop(0)

      # concatenation des positions du demi-cercle aux listes de positions x et y
      list_x += la_x
      list_y += la_y
      list_cap += lc

      # ajout des pas de temps du demi-cercle dans la liste de temps
      dts2 = [t + Temps[-1] for t in dts]
      Temps += dts2

      # calcul du cap après parcours du demi-cercle
      if orientation_virage != 'Gauche':
         new_cap = cap + 180
      else:
         new_cap = cap - 180

      # parcours de la seconde ligne droite
      x2,y2,rdts,rlc = rectiligne(list_x[-1],list_y[-1],v0,new_cap,tps_rect2-pas_de_temps,pas_de_temps)
      x2.pop(0)
      y2.pop(0)
      rdts.pop(0)
      rlc.pop(0)

      list_x += x2
      list_y += y2
      list_cap += rlc

      # ajout des pas de temps de la seconde ligne droite dans la liste de temps
      rdts2 = [t + Temps[-1] for t in rdts]
      Temps += rdts2

      return list_x , list_y , Temps , list_cap



def trajectoireU(x0,y0,v0,cap,rayon_courbure,angle_rotation,pas_de_temps,temps_acquisition,orientation_virage):
    """
    v0 = vitesse initiale en m/s
    rayon_courbure = rayon de courbure de l'arc du U, en m
    angle_rotation = angle de rotation de l'arc du U, en degrés
    pas_de_temps = pas de temps entre deux points successifs, en s
    temps_aquisition = temps total de l'acquisition, en s

   SORTIE : liste des positions en x
            liste des positions en y
            liste des temps
            liste des caps à tout instant
    """

    # initialisation des variables necessaires au parcours de l'arc
    R = rayon_courbure
    th = radians(angle_rotation)

    # calcul du temps necessaire au parcours de l'arc de cercle #
    tps_arc = th*R/v0

    if tps_arc > temps_acquisition:
        raise ValueError("Temps d'acquisition trop petit : impossible de faire l'arc de cercle voulu")
    
    elif tps_arc == temps_acquisition:
        return calculArc(x0,y0,v0,cap,R,pas_de_temps,angle_rotation,orientation_virage)
    
    else:
        # calcul du temps de parcours des lignes droites du U
        tps_restant = temps_acquisition - tps_arc
        tps_rect = tps_restant/2

        # parcours de la première ligne droite du U
        list_x , list_y , list_t , list_cap = rectiligne(x0,y0,v0,cap,tps_rect,pas_de_temps)

        # parcours de l'arc de cercle
        la_x , la_y , la1_t , lc = calculArc(list_x[-1],list_y[-1],v0,cap,R,pas_de_temps,angle_rotation,orientation_virage)
        la_x.pop(0)
        la_y.pop(0)
        la1_t.pop(0)
        lc.pop(0)

        # ajout des pas de temps de l'arc de cercle dans la liste de temps
        la2_t = [t + list_t[-1] for t in la1_t]        
        list_t += la2_t

        # concatenation de l'arc de cercle à la première ligne droite
        list_x += la_x
        list_y += la_y
        list_cap += lc

        # calcul du cap après parcours de l'arc de cercle
        if orientation_virage != 'Gauche':
            new_cap = cap + angle_rotation
        else:
            new_cap = cap - angle_rotation    


        # parcours de la seconde ligne droite du U
        lr2_x , lr2_y , lr2_t , lc2= rectiligne(la_x[-1],la_y[-1],v0,new_cap,tps_rect,pas_de_temps)
        lr2_x.pop(0)
        lr2_y.pop(0)
        lr2_t.pop(0)
        lc2.pop(0)

        # ajout des pas de temps de la seconde ligne droite dans la liste de temps
        lr2_2_t = [t + list_t[-1] for t in lr2_t]
        list_t += lr2_2_t

        # concatenation de l'arc de cercle à la première ligne droite
        list_x += lr2_x
        list_y += lr2_y
        list_cap += lc2

        return list_x , list_y , list_t , list_cap

###############
### Helpers ###
###############


def read_xml(data):
    
    # Parsage du fichier XML

    data = et.parse(data)
    root = data.getroot()

    # Récupération des données
    latitude = (root.find("position_initiale/latitude").text)
    longitude = (root.find("position_initiale/longitude").text)
    altitude = (root.find("altitude").text)
    cap = (root.find("cap").text)
    vitesse = (root.find("vitesse").text)
    pas_de_temps = (root.find("pas_de_temps").text)

    return latitude,longitude,altitude,vitesse,cap,pas_de_temps

def getXY(latitude,longitude):

    # Définition de la projection WGS84 = 'EPSG:4326'(système de coordonnées géodésiques)
    # Définition de la projection Mercator ='EPSG:3857'

    transformer = Transformer.from_crs('4326','3857')

    # Conversion des coordonnées cartésiennes en coordonnées géodésiques
    x0,y0 = transformer.transform(latitude,longitude)
    return x0,y0

def getLatLong(x,y):
    # Définition de la projection WGS84 = 'EPSG:4326'(système de coordonnées géodésiques)
    # Définition de la projection Mercator ='EPSG:3857'

    transformer = Transformer.from_crs('3857','4326') #

    # Conversion des coordonnées cartésiennes en coordonnées géodésiques
    latitude,longitude = transformer.transform(x,y)

    return latitude, longitude

def tracerTrajectoireCarte(X,Y,nom_carte):
    # Conversion des coordonnées (m) en latitude et longitude dans le système WGS84
    trajectoire=[[X[i],Y[i]] for i in range(len(X))]

    # Création de la carte centrée sur les coordonnées de départ
    carte = folium.Map(location=trajectoire[0], zoom_start=12)

    # Ajout de la trajectoire à la carte
    folium.PolyLine(locations=trajectoire, color='red', weight=2.5).add_to(carte)

    # Affichage de la carte
    carte.save(os.path.abspath(os.path.join(os.path.dirname(__file__), 'Cartes/'+nom_carte)))
    return 0

def sauvegarderResultats(X, Y, T, Z, cap, nom_fichier):
    dossier_resultats = "Résultats"  # Nom du dossier de résultats

    # Récupérer le chemin absolu du répertoire parent du fichier Python
    chemin_dossier_parent = os.path.dirname(os.path.abspath(__file__))

    # Créer le chemin complet du dossier de résultats
    chemin_dossier_resultats = os.path.join(chemin_dossier_parent, dossier_resultats)

    # Créer le dossier s'il n'existe pas déjà
    if not os.path.exists(chemin_dossier_resultats):
        os.makedirs(chemin_dossier_resultats)

    chemin_fichier = os.path.join(chemin_dossier_resultats, nom_fichier)

    with open(chemin_fichier, 'w') as fichier:
        fichier.write("#Fichier contenant sur chaque ligne : Latitude en degré, Longitude en degré, Altitude en mètre, Temps en secondes, Cap en degré\n")
        for i in range(len(X)):
            fichier.write((str(X[i]) +','+ str(Y[i]) +','+ str(Z[i]) +','+ str(T[i]) +','+ str(cap[i])).encode('ascii').decode('ascii') + "\n")

def telechargerResultats(X, Y, T, Z, cap, chemin_fichier):
    with open(chemin_fichier, 'w') as fichier:
        fichier.write("#Fichier contenant sur chaque ligne : Latitude en degré, Longitude en degré, Altitude en mètre, Temps en secondes, Cap en degré\n")
        for i in range(len(X)):
            fichier.write((str(X[i]) +','+ str(Y[i]) +','+ str(Z[i]) +','+ str(T[i]) +','+ str(cap[i])).encode('ascii').decode('ascii') + "\n")

def tracerCarteDynamique(X,Y,T,nom_carte):
    if len(X) != len(Y):
        print("ERREUR FONCTION tracerCarte")
        print("les vecteurs X et Y ne sont pas de la même taille")
        return 0
    
    else:
        # Conversion du temps
        T1 = []
        for utc in T:
            heures, minutes, secondes = convertir_secondes(utc) 
            dt = datetime(2023, 6, 14, heures, minutes, secondes)
            dt_str = dt.strftime("%Y-%m-%dT%H:%M:%S")
            T1.append(dt_str)

        # Conversion des coordonnées (m) en latitude et longitude dans le système WGS84
        trajectoire1 = [[Y[i], X[i]] for i in range(len(X))]

        # Création de la carte centrée sur les coordonnées de départ
        carte = folium.Map(location=[X[0], Y[0]], zoom_start=12)

        # Création de la ligne de trajectoire
        line = {
            "type": "Feature",
            "geometry": {
                "type": "LineString",
                "coordinates": trajectoire1,
            },
            "properties": {
                "times": T1,
                "style": {
                    "color": "red",
                    "weight": 2,
                },
            },
        }

        # Ajout de la ligne de trajectoire à la carte
        plugins.TimestampedGeoJson(
            {
                "type": "FeatureCollection",
                "features": [line],
            },
            period="PT1S",
            add_last_point=True,
            auto_play=True,  # Activer la lecture automatique
            loop=False,  # Désactiver la lecture en boucle
            max_speed=20,  # Vitesse maximale de 20 (peut être ajustée)
            min_speed=0.1,  # Vitesse minimale de 0.1 (peut être ajustée)
        ).add_to(carte)

        # Affichage de la carte
        carte.save(os.path.abspath(os.path.join(os.path.dirname(__file__), 'Cartes/'+nom_carte)))
        return 0
    
def convertir_secondes(secondes):
    heures = secondes // 3600
    reste = secondes % 3600
    minutes = reste // 60
    secondes = reste % 60

    return int(heures), int(minutes), int(secondes)    



###############################################
### Fonctions qui génèrent les trajectoires ###
###############################################


def generateRectiligne(lat0,long0,v0,cap,temps_acquisition,pas_de_temps):
    x0,y0 =getXY(lat0,long0)
    X,Y,T,liste_cap = rectiligne(x0,y0,v0,cap,temps_acquisition,pas_de_temps)
    Latitudes,Longitudes=getLatLong(X,Y)
    return Latitudes,Longitudes,T,liste_cap


def generateCirculaire(lat0,long0,v0,cap,rayon_courbure,temps_acquisition,pas_de_temps,orientation_virage):
    x0,y0 =getXY(lat0,long0)
    X,Y,T,liste_cap = circulaire(x0,y0,v0,cap,rayon_courbure,temps_acquisition,pas_de_temps,orientation_virage)
    Latitudes,Longitudes=getLatLong(X,Y)
    return Latitudes,Longitudes,T,liste_cap

def generateZigzag(lat0,long0,v0,cap,temps_acquisition,pas_de_temps,nbr_lignes,l_lignes,angle_chgt,virage_1,type_zigzag):
    x0,y0 =getXY(lat0,long0)
    X,Y,T,liste_cap = zigzag(x0,y0,v0,cap,temps_acquisition,pas_de_temps,nbr_lignes,l_lignes,angle_chgt,virage_1,type_zigzag)
    Latitudes,Longitudes=getLatLong(X,Y)
    return Latitudes,Longitudes,T,liste_cap

def generateZigzagAleatoire(lat0,long0,v0,cap,temps_acquisition,pas_de_temps):
    x0,y0 =getXY(lat0,long0)
    X,Y,T,liste_cap = zigzag_full_alea(x0,y0,v0,cap,temps_acquisition,pas_de_temps)
    Latitudes,Longitudes=getLatLong(X,Y)
    return Latitudes,Longitudes,T,liste_cap

def generateHippodrome(lat0,long0,v0,cap,rayon_courbure,temps_acquisition,pas_de_temps,orientation_virage):
    x0,y0 = getXY(lat0,long0)
    X,Y,T,liste_cap = hippodrome(x0,y0,v0,cap,rayon_courbure,temps_acquisition,pas_de_temps,orientation_virage)
    Latitudes,Longitudes=getLatLong(X,Y)
    return Latitudes,Longitudes,T,liste_cap

def generateArc(lat0,long0,v0,cap,rayon_courbure,temps_acquisition,pas_de_temps,angle_rot,orientation_virage):
    x0,y0 = getXY(lat0,long0)
    X,Y,T,liste_cap = trajectoireArc(x0,y0,v0,cap,rayon_courbure,temps_acquisition,pas_de_temps,angle_rot,orientation_virage)
    Latitudes,Longitudes=getLatLong(X,Y)
    return Latitudes,Longitudes,T,liste_cap

def generateU(lat0,long0,v0,cap,rayon_courbure,angle_rot,pas_de_temps,temps_acquisition,orientation_virage):
    x0,y0 = getXY(lat0,long0)
    X,Y,T,liste_cap = trajectoireU(x0,y0,v0,cap,rayon_courbure,angle_rot,pas_de_temps,temps_acquisition,orientation_virage)
    Latitudes,Longitudes=getLatLong(X,Y)
    return Latitudes,Longitudes,T,liste_cap




