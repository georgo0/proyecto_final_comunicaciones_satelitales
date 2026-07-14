import math
from datetime import datetime, timedelta
from zoneinfo import ZoneInfo

#Lectura de datos 
def parser(ruta): #lee archivo txt que contiene coordenadas y TLE
    try:
        with open(ruta, 'r') as file:
            lines = file.readlines()      #Leer archivo linea por linea
            line1 = lines[0]              #Obtener primera linea de Tle
            line2 = lines[1]              #Obtener segunda linea de Tle
            
            #se extrae el dato de la posiciones de columna especifica
            #indicado en el formato estandar de TLE
            data = {
                "catalog_number": line1[2:7].strip(),
                "classification": line1[7:8].strip(),
                "epoch_year": line1[18:20],
                "epoch_day": line1[20:32],
                "ascencion_recta_nodo": line2[18:25].strip(),
                "arg_de_perigeo": line2[34:42].strip(),
                "inclinacion": line2[8:16].strip(),
                "excentricidad": f"0.{line2[26:33]}".strip(),
                "anomalia_media": line2[43:51].strip(), #en grados
                "medium_move": line2[52:69].strip(),  # en rev/por dia
                "obs_lattle" : lines[2].strip(), #latitud observador
                "obs_long" : lines[3].strip() #longitud observador
            }  
            return data

    except FileNotFoundError:
        print(f"Error: The file '{ruta}' was not found.")
        
#Obtener Elementos Keplerianos ------------------------------------------------------------------------------------      
def get_anomalia_excentrica(anommedia,ex,max=15): #anomalia excentrica en grados
    M = math.radians(float(anommedia))            #Conversion de grados a radianes
    exc = float(ex)                               #Ingresar excentricidad de la orbita
    tolerance = math.pow(10,-8)                   #Tolerancia como 10 a la -8 
    Aexc = M                                      #Valor inicial a iterar
    for i in range(max):
        #Metodo iterativo Newton-Raphson
        f = Aexc - exc * math.sin(Aexc) - M    #Funcion de kepler
        df = 1 - exc * math.cos(Aexc)          #Derivada de funcion de kepler 
        Aexc_next = Aexc - (f/df)              #Calculo de siguiente iteracion 
        
        if abs(Aexc-Aexc_next) < tolerance:    #Si la iteracion actual esta bajo tolerancia
            print("Aexc bajo tolerancia")      
            return math.degrees(Aexc_next)    
        
        Aexc = Aexc_next                      #se actualiza el valor de la iteracion
        
    print("Numero de Iteraciones excedidas")
    return math.degrees(Aexc)                 

def get_anomalia_verdadera(anom_exc,ex):            #calcula la posicion real en la elipse
    Aexc = math.radians(float(anom_exc))            #pasar a radianes para las funciones trigonométricas
    exc = float(ex)                                 #asegurar que la excentricidad sea float
    
    num = math.sqrt(1 + exc) * math.sin(Aexc / 2)   #parte de arriba de la formula
    den = math.sqrt(1 - exc) * math.cos(Aexc / 2)   #parte de abajo de la formula
    
    halfv = math.atan2(num, den)                    #arcotangente para sacar la mitad del angulo (v/2)
    anom_true = 2 * halfv                           #multiplicar por 2 para tener la anomalia completa

    return math.degrees(anom_true)                  #se devuelve en grados para los prints

def get_axis_a(medium_move_dias): 
    mov_medio = float(medium_move_dias)*(2*math.pi/86400) #pasar de revoluciones/dia a radianes/segundo
    u = 398600.4418                                       #constante gravitacional de la tierra en km
    axisa = (u / (mov_medio**2))**(1/3)                   #formula para despejar el semieje mayor (a)
    return axisa

def distancia_radial(a,e,anom_true): #calcula la distancia desde el centro de la tierra al satelite
    axis = a
    exc = float(e)
    Atrue = math.radians(anom_true)      #necesitamos la anomalia verdadera en radianes
    
    num = axis * (1 - exc**2)            #numerador de la ecuacion de trayectoria polar
    denom = 1 + exc * math.cos(Atrue)    #denominador
    radial_distance = num/denom          #distancia en km
    
    return radial_distance
    
#Obtencion de angulo GMST --------------------------------------------------------------------------------------  
def get_time(): #obtiene la hora actual del computador
    zona_local = ZoneInfo("America/Santiago") #fijar la zona horaria a Chile
    hora_local = datetime.now(zona_local)     #sacar la hora exacta de ahora
    hora_utc = hora_local.astimezone(ZoneInfo("UTC")) #convertirla al estandar UTC
    return hora_utc
    
def get_julian_datetime(date): #convierte fecha normal a dias julianos (formato continuo)
    #En Enero o Febrero se traspasa un Año a 12 meses para que la formula funcione
    if date.month <= 2:
        date.year -= 1
        date.month += 12
        
    #Formula para sacar los dias julianos basada en el año, mes, dia y hora
    julian_datetime = 367 * date.year - int((7 * (date.year + int((date.month + 9) / 12.0))) / 4.0) \
    + int((275 * date.month) / 9.0) + date.day + 1721013.5 + (date.hour + date.minute / 60.0 + \
    date.second / math.pow(60,2)) / 24.0
    
    return julian_datetime
    
def gmst_from_jd(jd): #saca el angulo de rotacion de la tierra en Greenwich
    t = (jd - 2451545.0) / 36525.0 #siglos julianos transcurridos desde el año 2000
    
    # Ecuacion polinomial obtenida de pag 215 de Libro Fundamentals of Astrodynamics
    theta_gmst_sec = (
        67310.54841
        + (876600.0 * 3600.0 + 8640184.812866) * t
        + 0.093104 * t**2
        - 6.2e-6 * t**3
    )

    # El modulo sirve para que si la tierra dio más de una vuelta, el angulo vuelva a estar entre 0 y 86400 seg
    theta_gmst_sec = math.fmod(theta_gmst_sec, 86400.0)
    if theta_gmst_sec < 0:
        theta_gmst_sec += 86400.0

    gmst_deg = theta_gmst_sec / 240.0  #dividir para pasar de segundos de tiempo a grados de rotacion
    return gmst_deg
    
def deg_to_hms(deg): #funcion extra por si queremos ver el angulo en formato horas, minutos, segundos
    total_hours = deg / 15.0 # 15 grados es 1 hora de rotacion
    hours = int(total_hours)
    minutes = int((total_hours - hours) * 60)
    seconds = (total_hours - hours - minutes/60) * 3600
    return hours, minutes, seconds

def gmst_hms_from_jd(jd): #conecta las dos funciones anteriores
    gmst_deg = gmst_from_jd(jd)
    return deg_to_hms(gmst_deg)

#Coordenadas ECI -------------------------------------------------------------------------------
def get_eci(inc,arg_perig,anom_verdadera,asc_node,ra): #ubica al satelite en el plano 3D ECI
  #se pasan todos los angulos a radianes si o si para hacer math.cos y math.sin
  inclinacion_rad = math.radians(float(inc))
  asc_node_rad = math.radians(float(asc_node))
  arg_perig_rad = math.radians(float(arg_perig))
  anom_verd_rad = math.radians(float(anom_verdadera))
  distancia = float(ra) #distancia radial ya sacada en funcion anterior
    
  #Coordenadas en el plano orbital (2D, asumiendo que la orbita esta plana)
  x_pqw = distancia * math.cos(anom_verd_rad)
  y_pqw = distancia * math.sin(anom_verd_rad)

  #Precalcular los senos y cosenos para no recargar la matriz abajo
  sin_i = math.sin(inclinacion_rad)
  cos_i = math.cos(inclinacion_rad)
  sin_O = math.sin(asc_node_rad)
  cos_O = math.cos(asc_node_rad)
  sin_w = math.sin(arg_perig_rad)
  cos_w = math.cos(arg_perig_rad)

  #Matriz de rotacion (Toma el plano 2D y lo rota segun la inclinacion y los nodos)
  x_eci = x_pqw * (cos_O * cos_w - sin_O * sin_w * cos_i) - y_pqw * (cos_O * sin_w + sin_O * cos_w * cos_i)
  y_eci = x_pqw * (sin_O * cos_w + cos_O * sin_w * cos_i) - y_pqw * (sin_O * sin_w - cos_O * cos_w * cos_i)
  z_eci = x_pqw * (sin_w * sin_i) + y_pqw * (cos_w * sin_i)
  
  sat_vector = [x_eci,y_eci,z_eci] #Vector final X Y Z del satelite
  return sat_vector
  
def posicion_obs(longitud, latitud, gmst_deg): #ubica a nuestra antena en el plano 3D ECI
  # Constantes WGS84 para no asumir que la tierra es una esfera perfecta (es achatada)
  radio_ecuatorial = 6378.137  # en kilometros
  achatamiento = 1 / 298.257223563 
  e2 = 2 * achatamiento - achatamiento**2  # Cuadrado de la excentricidad del planeta

  #pasar a radianes
  lat_rad = math.radians(float(latitud))
  lon_rad = math.radians(float(longitud))
  gmst_rad = math.radians(float(gmst_deg))

  # Radio de curvatura vertical en nuestra latitud
  N = radio_ecuatorial / math.sqrt(1 - e2 * math.sin(lat_rad)**2)

  # Posicion ECEF (Sistema de coordenadas que gira junto con la tierra)
  X_ecef = N * math.cos(lat_rad) * math.cos(lon_rad)
  Y_ecef = N * math.cos(lat_rad) * math.sin(lon_rad)
  Z_ecef = (N * (1 - e2)) * math.sin(lat_rad)

  # Rotar ECEF a ECI usando el angulo GMST (segun la hora del dia)
  x_obs_eci = X_ecef * math.cos(gmst_rad) - Y_ecef * math.sin(gmst_rad)
  y_obs_eci = X_ecef * math.sin(gmst_rad) + Y_ecef * math.cos(gmst_rad)
  z_obs_eci = Z_ecef

  obs_vector = [x_obs_eci,y_obs_eci,z_obs_eci] #Vector X Y Z de nuestra casa
  return obs_vector

def get_view_vector(sat_vector, obs_vector): #Linea de vision desde la antena al satelite
  x = sat_vector[0] - obs_vector[0] #Resta simple coordenada a coordenada
  y = sat_vector[1] - obs_vector[1]
  z = sat_vector[2] - obs_vector[2]
  view_vector = [x,y,z] #El nuevo origen es la antena
  return view_vector
  
#Coordenadas Horizontales ----------------------------------------------------------------------------------

def transf_a_sez(latitud, longitud, angulo_horario, r_rel): #Transforma el vector ECI al horizonte local (Sur Este Cenit)
    r_eci = r_rel
    lat_deg = float(latitud)
    lon_deg = float(longitud)
    gst_deg = float(angulo_horario) #GMST
    
    # Tiempo Sidereo Local (LST) = Greenwich + Longitud local
    lst_deg = gst_deg + lon_deg
    
    # angulos para hacer la rotacion
    phi = math.radians(lat_deg)
    theta = math.radians(lst_deg)
    
    # desempaquetar el vector de vision
    dx = r_eci[0]
    dy = r_eci[1]
    dz = r_eci[2]
    
    # Multiplicacion de la matriz de rotacion SEZ (topocentrica)
    S = dx * math.sin(phi) * math.cos(theta) + dy * math.sin(phi) * math.sin(theta) - dz * math.cos(phi)
    E = -dx * math.sin(theta) + dy * math.cos(theta)
    Z = dx * math.cos(phi) * math.cos(theta) + dy * math.cos(phi) * math.sin(theta) + dz * math.sin(phi)
    
    return [S, E, Z]

def calc_azimut_elev(vector): #Convierte las coordenadas XYZ locales a angulos para apuntar
    S = vector[0] #Coordenada Sur
    E = vector[1] #Coordenada Este
    Z = vector[2] #Coordenada Cenit (arriba)
    
    # Elevacion: arco tangente del eje Z divido por la hipotenusa de S y E
    elevacion_rad = math.atan2(Z, math.sqrt(S**2 + E**2))
    elevacion_deg = math.degrees(elevacion_rad)
    
    # Azimut: arco tangente del Este contra el Sur negativo (Norte)
    azimuth_rad = math.atan2(E, -S)
    azimuth_deg = math.degrees(azimuth_rad)
    
    # El azimut debe ser un compas de 0 a 360, si da negativo se suma 360
    if azimuth_deg < 0:
        azimuth_deg += 360

    return azimuth_deg, elevacion_deg

#Ejecucion Secuencial-----------------------------------------------------------------------------------------
# Cambiar archivo de input, usado como referencia NOAA21
ruta = 'noatle.txt' 
info = parser(ruta)

horautc = get_time() #hora para iniciar todos los calculos
dias_julianos = get_julian_datetime(horautc)
ang_gmst = gmst_from_jd(dias_julianos) #angulo de la tierra
hora_gmst = gmst_hms_from_jd(dias_julianos)

anio_str = info["epoch_year"]
#Si el año es menor a 57 es del 2000, si no es de 1900 (logica de NORAD)
anio = 2000 + int(anio_str) if int(anio_str) < 57 else 1900 + int(anio_str)
dias_fraccionales = float(info["epoch_day"])

#Se crea un objeto fecha con la epoca original del TLE
fecha_base = datetime(anio, 1, 1, tzinfo=ZoneInfo("UTC"))
fecha_utc = fecha_base + timedelta(days=dias_fraccionales - 1)

#Delta T: cuantos dias exactos pasaron entre el TLE y la hora de recien
delta_t_dias = (horautc - fecha_utc).total_seconds() / 86400.0

#Se actualiza la anomalia media multiplicando los dias por el movimiento diario
M_tle = float(info["anomalia_media"])
n_grados_dia = float(info["medium_move"]) * 360.0
M_actual_grados = (M_tle + (n_grados_dia * delta_t_dias)) % 360.0

#Cascada de calculos keplerianos
Aexcentrica = get_anomalia_excentrica(M_actual_grados,info["excentricidad"])
Averdadera = get_anomalia_verdadera(Aexcentrica,info["excentricidad"])
axisa = get_axis_a(info["medium_move"])
dist = distancia_radial(axisa,info["excentricidad"],Averdadera)

#Calculo de vectores de posicion ECI
vector_sat = get_eci(info["inclinacion"],info["arg_de_perigeo"],Averdadera,info["ascencion_recta_nodo"],dist)
vector_obs = posicion_obs(info["obs_lattle"],info["obs_long"],ang_gmst)
vector_vista = get_view_vector(vector_sat,vector_obs)

#Paso final a horizonte local
vector_sez = transf_a_sez(info["obs_lattle"],info["obs_long"],ang_gmst,vector_vista)
vector_horiz = calc_azimut_elev(vector_sez)

#Muestra de resultados por pantalla
print("Resultados de Posicionamiento Satelital:")
print("---------------------")
print(f"Hora de calculo UTC: {horautc}")
print(f"Anomalia Media leida del TLE: {info['anomalia_media']}")
print(f"Anomalia Media propagada al dia de hoy: {M_actual_grados:.4f}")
print("---------------------")
print(f"Distancia radial del satelite: {dist:.2f} km")
print(f"Vector Observador ECI: {vector_obs}")
print(f"Vector Satelite ECI: {vector_sat}")
print("---------------------")
print(f"Elevacion final: {vector_horiz[1]:.2f} grados")
print(f"Azimut final: {vector_horiz[0]:.2f} grados")