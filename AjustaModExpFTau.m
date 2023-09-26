%*******************************************************************************
% Función: [s0,c,Exx,du] = AjustaModExpFTau(datos,tau,saltos,pesos,eg)
%          [s0,c,Exx,du] = AjustaModExpFTau(datos,tau,saltos,pesos)
%          [s0,c,Exx,du] = AjustaModExpFTau(datos,tau,saltos)
%          [s0,c,Exx,du] = AjustaModExpFTau(datos,tau)
%
% Propósito: Ajusta una serie de valores de deformación al modelo exponencial
%            s=s0*(1-exp(-t/tau))+Sum[c(t>=t_k)], donde los parámetros a
%            determinar son s0 y c(t>=t_k), generando estos últimos una función
%            a trozos donde s0 es comun.
%
% Entrada: - datos: Matriz de dos o tres columnas con los datos de trabajo:
%                   - Col. 1: Instante de tiempo.
%                   - Col. 2: Valor observado de s.
%                   - Col. 3: Desviación típica del valor observado (esta
%                             columna puede estar presente o no).
%          - tau: Valor del parámetro tau.
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
%
% Salida: - s0: Valor ajustado del parámetro s0.
%         - c: Matriz de dos filas:
%              - Fil. 1: Instantes de tiempo de los saltos.
%              - Fil. 2: Valores ajustados de los posibles saltos, en el mismo
%                        orden que hayan sido pasados en el argumento 'saltos'.
%         - Exx: Matriz de varianzas-covarianzas a posteriori de las incógnitas,
%                en el orden s0, tau y c(2,:).
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
% Historia: 15-09-2023: Creación de la función
%           José Luis García Pallero, jgpallero@gmail.com
%*******************************************************************************
function [s0,c,Exx,du] = AjustaModExpFTau(datos,tau,saltos,pesos,eg)

%Posibles valores por omisión
if nargin<5
    %Por omisión no se realiza detección de errores groseros
    eg = 0;
    if nargin<4
        %Trabajo sin ponderar
        pesos = 0;
        if nargin<3
            %No considero saltos
            saltos = [];
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
    %Me quedo sólo con los saltos posteriores al día de inicio
    pos = saltos>=min(datos(:,1));
    saltos = saltos(pos);
    %Número real de saltos y número máximo y mínimo de días
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
%En principio se usan todos los puntos
du = logical(ones(fil,1));
%Entro en un bucle infinito
while 1
    %Realizo el ajuste
    [s0,c,v,A,P,Qxx] = AjustaModExpFTauAux(datos,tau,c,pesos);
    %Grados de libertad del ajuste
    gl = fil-(1+size(c,2));
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
%Función auxiliar para realizar el ajuste
function [s0,c,v,A,P,Q] = AjustaModExpFTauAux(datos,tau,c,pesos)

%Número de datos
fil = size(datos,1);
col = size(datos,2);
%Número de saltos
ns = size(c,2);
%Matrices de diseño y de pesos
A = zeros(fil,1+ns);
if pesos&&(col>2)
    P = diag(1.0./(datos(:,3).^2));
else
    P = eye(fil);
end
%Vector término independiente para el ajuste
l = datos(:,2);
%Primera columna de la matriz de diseño
A(:,1) = 1.0-exp(-datos(:,1)./tau);
%Recorro los posibles saltos
for i=1:ns
    %Posiciones de los puntos tras el salto
    pos = datos(:,1)>=c(1,i);
    %Columnas de la matriz de diseño correspondientes a los saltos
    A(pos,i+1) = 1.0;
end
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
%Residuos
v = A*x-l;
%Incógnitas
s0 = x(1);
if ns
    c(2,:) = x(2:end)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (c) 2023, J.L.G. Pallero. All rights reserved.
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
