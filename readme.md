Taller sobre *Implementación de procesos puntuales espacio-temporales con NIMBLE* impartido en el *IV Congreso & XV JORNADAS de Usuarios de R*.

La carpeta principal contiene los códigos para ejecutar los diferentes modelos descritos durante el taller, incluyendo la carga y preparación de datos, la definición de constantes y valores iniciales, y el ajuste de modelos con nimble. En concreto:

- "0 Data exploration.R" contiene algunas líneas de código para una exploración inicial de los datos.

- "1 Homogeneous Poisson model.R" contiene el código para ajustar el modelo de Poisson homogéneo.

- "2 Inhomogeneous Poisson model v1.R", "2 Inhomogeneous Poisson model v2.R", "2 Inhomogeneous Poisson model v3.R" y "2 Inhomogeneous Poisson model v4.R" contienen el código para ajustar las cuatro versiones siguientes de un modelo de Poisson no homogéneo: incluyendo covariables espaciales (v1), incluyendo covariables espaciales y un efecto temporal que decae a partir del día t = 200 (v2), incluyendo covariables espaciales y un random walk de primer orden para el efecto temporal (v3), e incluyendo covariables espaciales y un término periódico (trigonométrico) para el efecto temporal (v4).

- "3 Splines for the spatial intensity.R" contiene el código para reemplazar la estimación de la intensidad espacial mediante covariables por una estimación basada en funciones spline definidas sobre la ventana espacial.

- "4 Residual analysis.R" contiene el código para realizar un análisis residual espacial de cualquiera de los modelos ajustados (ten en cuenta que el resultado del modelo debe estar en la carpeta "Outputs", por lo que debes ejecutar el código "1 Modelo Poisson homogéneo.R" o cualquier otro de los códigos para realizar un análisis residual).

Además:

- La carpeta **Models** contiene los códigos nimble de cada uno de los modelos bayesianos de procesos puntuales explicados.

- La carpeta **Data** contiene los datos necesarios para utilizar los códigos. Ten en cuenta que, para evitar violar la privacidad de los datos, las ubicaciones temporales se han obtenido a partir de una distribución uniforme en [0,365] (los resultados que se muestran en las diapositivas corresponden al conjunto de datos real).

- La carpeta **Outputs** debe utilizarse para guardar los resultados de los modelos ajustados según los códigos. La carpeta aparece vacía porque los resultados son bastante grandes. Ten cuidado con el número de iteraciones/cadenas que eliges para ajustar el modelo.

- La carpeta **Slides** contiene las diapositivas empleadas durante el taller.
