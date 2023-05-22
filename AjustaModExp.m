%*******************************************************************************
% Función: [s0,tau,c,Exx,iter,du] = AjustaModExp(datos,saltos,pesos,eg,parada)
%          [s0,tau,c,Exx,iter,du] = AjustaModExp(datos,saltos,pesos,eg)
%          [s0,tau,c,Exx,iter,du] = AjustaModExp(datos,saltos,pesos)
%          [s0,tau,c,Exx,iter,du] = AjustaModExp(datos,saltos)
%          [s0,tau,c,Exx,iter,du] = AjustaModExp(datos)
%
% Propósito: Ajusta una serie de valores de deformación al modelo exponencial
%            s=s0*(1-exp(-t/tau))+Sum[c(t>=t_k)], donde los parámetros a
%            determinar son s0, tau, y c(t>=t_k), generando estos últimos una
%            función a trozos donde s0 y tau son comunes a todos.
%
% Entrada: - datos: Matriz de dos o tres columnas con los datos de trabajo:
%                   - Col. 1: Instante de tiempo.
%                   - Col. 2: Valor observado de s.
%                   - Col. 3: Desviación típica del valor observado (esta
%                             columna puede estar presente o no).
%          - saltos: Vector con los instantes de tiempo t_k en que se considera
%                    que ha habido un salto en el valor de s para t>=t_k. Puede
%                    pasarse un vector vacío (que es el valor por omisión).
%          - pesos: Escalar para indicar si se realiza un ajuste ponderado o no:
%                   - 0: El ajuste se realiza sin ponderar (valor por omisión)
%                   - Distinto de 0: El ajuste se realiza ponderado, siendo los
%                                    pesos iguales al inverso de las varianzas
%                                    de cada observación (si no están presentes
%                                    en 'datos' no se pondera).
%          - eg: Identificador para la detección de errores groseros:
%                - 0 o >= 1: No se realiza ningún test de detección de errores.
%                - Menor que 0: Se eliminan los puntos (uno a uno) cuyo residuo
%                               sea mayor que el valor pasado.
%                - Número en el intervalo (0, 1): Indica el nivel de
%                  significación para realzar el test de Pope, donde se eliminan
%                  los puntos de uno en uno.
%                Por omisión este parámetro vale 0.
%          - parada: Vector de dos elementos con los criterios de parada del
%                    proceso iterativo de ajuste:
%                    - Pos. 1: Criterio de parada en función de la relación
%                              norm(dx)/norm(x), donde norm() es la norma L2, dx
%                              es la corrección al vector de incógnitas de la
%                              iteración actual y x el vector de incógnitas.
%                    - Pos. 2: Número máximo de iteraciones.
%                    Por omisión, este argumento vale [0.0001 10].
%
% Salida: - s0: Valor ajustado del parámetro s0.
%         - tau: Valor ajustado del parámetro tau.
%         - c: Matriz de dos filas:
%              - Fil. 1: Instantes de tiempo de los saltos.
%              - Fil. 2: Valores ajustados de los posibles saltos, en el mismo
%                        orden que hayan sido pasados en el argumento 'saltos'.
%         - Exx: Matriz de varianzas-covarianzas a posteriori de las incógnitas,
%                en el orden s0, tau y c(2,:).
%         - iter: Número de iteraciones realizadas.
%         - du: Vector booleano de la misma longitud que las filas de 'datos'.
%               Los posibles valores de cada elemento son:
%               - 0: El elemento correspondiente de la fila de 'datos' no se ha
%                    utilizado en el ajuste por considerarse error grosero.
%               - 1: El elemento correspondiente de la fila de 'datos' se ha
%                    utilizado en el ajuste.
%         - s02: Varianza de referencia del observable de peso unidad a
%                posteriori para el ajuste (s02=v'*P*v/gl, donde gl son los
%                grados de libertad, filas-columnas de la matriz de diseño).
%
% Historia: 23-02-2023: Creación de la función
%           J.L.G. Pallero and M. Charco, jgpallero@gmail.com, m.charco@csic.es
%           17-04-2023: Adición del parámetro de salida s02
%           J.L.G. Pallero, jgpallero@gmail.com
%*******************************************************************************
function [s0,tau,c,Exx,iter,du,s02] = AjustaModExp(datos,saltos,pesos,eg,parada)

%Posibles valores por omisión
if nargin<5
    %Impongo el criterio de parada
    parada = [0.0001 10];
    if nargin<4
        %Por omisión no se realiza detección de errores groseros
        eg = 0;
        if nargin<3
            %Trabajo sin ponderar
            pesos = 0;
            if nargin<2
                %No considero saltos
                saltos = [];
            end
        end
    end
end
%Copia de la matriz de datos originales
datos_orig = datos;
%Compruebo si no hay saltos e inicializo la variable de salida correspondiente
if length(saltos)==0
    %Número de saltos
    ns = 0;
    c = zeros(2,0);
else
    %Número real de saltos y número máximo de días
    ns = 0;
    nmd = max(datos(:,1));
    %Recorro el número de saltos
    for i=1:length(saltos)
        %Compruebo si hay datos para llegar al salto
        if nmd>=saltos(i)
            %Aumento el número de saltos
            ns = ns+1;
        end
    end
    %Dimensiono la matriz de saltos
    c = zeros(2,ns);
    c(1,:) = saltos(1:ns);
end
%Dimensiones de la matriz de datos
fil = size(datos,1);
%Posiciones de los 1/6 últimos datos y de los 1/15 primeros
nu = floor(fil/6.0);
pu = (fil-nu):fil;
np = ceil(fil/15.0);
pp = 1:np;
%Primera aproximación de s0 como la media de los 1/6 últimos datos
s0 = mean(datos(pu,2));
%La primera aproximación a tau la calculo con el valor de s0 inicial y el
%acortamiento medio del instante medio de los primeros 1/15 puntos
%Hay que tener cuidado para que los valores de tiempo y acortamiento no sean
%cero, lo que dará un valor inicial de rtau igual a NaN
while 1
    %Valor inicial de tau
    tau = mean(datos(pp,1))/log(s0/(s0-mean(datos(pp,2))));
    %Comprobamos si el valor obtenido es NaN y si se pueden seguir añadiendo
    %puntos a los primeros
    if isnan(tau)&&(length(pp)<fil)
        %Añado un punto más a los primeros para estimar el valor inicial de tau
        pp = [pp pp(end)+1];
    else
        %Salimos del bucle
        break;
    end
end
%En principio se usan todos los puntos
du = logical(ones(fil,1));
%Entro en un bucle infinito
while 1
    %Realizo el ajuste
    [s0,tau,c,v,A,P,Qxx,iter] = AjustaModExpAux(datos,s0,tau,c,pesos,parada);
    %Grados de libertad del ajuste
    gl = fil-(2+size(c,2));
    %Varianza del observable de peso unidad a posteriori
    s02 = v'*P*v/gl;
    %Comprobamos si hay que realizar test de errores groseros
    if (eg<0.0)||((eg>0.0)&&(eg<1.0))
        if (eg>0.0)&&(eg<1.0)
            %Diagonal de A*Qxx*A'
            aqat = diag(A*Qxx*A');
            %Matriz cofactor a priori
            Q = 1.0./diag(P);
            %Matriz para almacenar el estadístico del test de Pope
            T = zeros(fil,1);
            %Calculamos el estadístico para el test
            T(:,1) = abs(v)./(sqrt(s02)*sqrt(Q-aqat));
            %Calculamos el nivel de significación para tau
            alfa0 = 1.0-(1.0-eg)^(1.0/fil);
            %Calculamos el valor de la T de Student correspondiente
            t = abs(tinv(alfa0/2.0,gl-1));
            %Calculamos el valor correspondiente de tau de Pope
            tauP = sqrt(gl)*t/sqrt(gl-1.0+t^2);
            %Puntos cuyo residuo supera el límite
            fuera = abs(v)>tauP;
        else
            %Buscamos los residuos mayores que el corte
            fuera = abs(v)>abs(eg);
        end
        %Compruebi si hay algún punto a eliminar
        if sum(fuera)
            %Busco los posibles puntos a eliminar
            datos_fuera = [datos abs(v)];
            datos_fuera = datos_fuera(fuera,:);
            [vmax,imax] = sort(datos_fuera(:,end),'descend');
            idmax = datos_fuera(imax(1),1);
            %Busco la posición del identificador en la matriz original
            pos = datos_orig(:,1)==idmax;
            %Quito el punto
            du(pos) = logical(0);
            datos = datos_orig(du,:);
            %Nuevo número de datos
            fil = size(datos,1);
            %Compruebo si hay que suprimir algún salto
            for i=size(c,2):-1:1
                pos = datos(:,1)>c(1,i);
                if sum(pos)==0
                    %Borro la columna
                    c(:,i) = [];
                end
            end
        else
            %Salimos del bucle
            break;
        end
    else
        %Salimos del bucle
        break;
    end
end
%Matriz de varianzas-covarianzas a posteriori
Exx = s02*Qxx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Función auxiliar para realizar el ajuste (ha de recibir los valores iniciales
%de las incógnitas)
function [s0,tau,c,v,A,P,Q,iter] = AjustaModExpAux(datos,s0,tau,c,pesos,parada)

%Número de datos
fil = size(datos,1);
col = size(datos,2);
%Número de saltos
ns = size(c,2);
%Matrices de diseño y de pesos
A = zeros(fil,2+ns);
if pesos&&(col>2)
    P = diag(1.0./(datos(:,3).^2));
else
    P = eye(fil);
end
%Vector de parámetros para el criterio de parada
x0 = [s0 tau c(2,:)]';
%Entramos en un bucle con, en principio, el número máximo de vueltas
for iter=1:parada(2)
    %Vector para almacenar la señal con los parámetros actuales
    sc = zeros(fil,1);
    %Calculo la señal en principio para todos los puntos
    sc = s0*(1.0-exp(-datos(:,1)./tau));
    %Recorro los posibles saltos
    for i=1:ns
        %Posiciones de los puntos tras el salto
        pos = datos(:,1)>=c(1,i);
        %Señal de los puntos posteriores al salto
        sc(pos) = sc(pos)+c(2,i);
        %Columnas de la matriz de diseño correspondientes a los saltos
        A(pos,i+2) = 1.0;
    end
    %Vector término independiente para el ajuste
    l = datos(:,2)-sc;
    %Primera columna de la matriz de diseño (derivada con respecto a s0)
    A(:,1) = 1.0-exp(-datos(:,1)./tau);
    %Segunda columna de la matriz de diseño (derivada con respecto a tau)
    A(:,2) = -s0*datos(:,1)./(tau^2).*exp(-datos(:,1)./tau);
    %Sistema normal
    N = A'*P*A;
    d = A'*P*l;
    %Solución
    if exist('OCTAVE_VERSION')
        Q = cholinv(N);
    else
        Q = inv(N);
    end
    x = Q*d;
    %Actualizo los parámetros
    s0 = s0+x(1);
    tau = tau+x(2);
    if ns
        c(2,:) = c(2,:)+x(3:end)';
    end
    %Actualizo el vector de parámetros para el criterio de parada
    x0 = [s0 tau c(2,:)]';
    %Criterio de parada
    if (norm(x)/norm(x0))<parada(1)
        %Salgo del bucle
        break;
    end
end
%Residuos
v = A*x-l;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (c) 2023, J.L.G. Pallero and M. Charco. All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%- Redistributions of source code must retain the above copyright notice, this
%  list of conditions and the following disclaimer.
%- Redistributions in binary form must reproduce the above copyright notice,
%  this list of conditions and the following disclaimer in the documentation
%  and/or other materials provided with the distribution.
%- Neither the name of the copyright holders nor the names of its contributors
%  may be used to endorse or promote products derived from this software without
%  specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
%INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
%OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
%ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
