%Carga de datos
datos = load('acortamiento-mazo-lp01.txt');
%Días a usar en los ajustes
dini = 15;
dfin = 125;
%Incrementos de presiones
dp = [0.01 0.02 0.025 0.05 0.1 0.2];
%Saltos, pesos y criterios de parada
% saltos = [0 45 49 62];
saltos = [62];
% saltos = [];
pesos = 1;
parada = [0.0001 10];
%Detección de errores groseros (eg==0 o eg>=1 no los busca, valor en (0,1)
%quiere decir el nivel de significación del test de Pope y un número negativo es
%el residuo máximo)
% eg = 0.01;
eg = -5.0;
% eg = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Posiciones de los días a usar y para dibujar
p_dini = find(datos(:,1)==dini);
if length(p_dini)==0
    error('El día indicado en la variable ''dini'' no existe en el fichero');
end
p_dfin = find(datos(:,1)==dfin);
if length(p_dfin)==0
    error('El día indicado en la variable ''dfin'' no existe en el fichero');
end
if p_dfin<p_dini
    error(['El valor de la variable ''dfin'' ha de ser mayor que al de la ',...
           'variable''dini''']);
end
%Duraciones estimadas y sus incertidumbres
d = zeros(length(dp),length(p_dini:p_dfin));
dt = zeros(length(dp),length(p_dini:p_dfin));
%Recorro los incrementos de presión
leyenda = {};
for i=1:length(dp)
    %Índice para las columnas de la matriz de almacenamiento de datos
    pos = 1;
    %Recorro las duraciones
    for j=p_dini:p_dfin
        %Realizo el ajuste
        [s0,tau,c,Exx,iter] = AjustaModExp(datos(1:j,:),saltos,pesos,eg,parada);
        if iter==parada(2)
            warning(['Se ha llegado al número máximo de iteraciones para ',...
                     'DeltaP=%f y %d datos . Se omite este ajuste en el ',...
                     'dibujo'],dp(i),length(1:j));
            %Almaceno los valores NaN
            d(i,pos) = NaN;
            dt(i,pos) = NaN;
        else
            %Calculo la duración y su incertidumbre
            ldp = log(dp(i));
            t = -tau*ldp;
            stdT = sqrt(ldp^2*Exx(2,2));
            %La almaceno
            d(i,pos) = t;
            dt(i,pos) = stdT;
        end
        %Actualizo el contador de posiciones
        pos = pos+1;
    end
    %Figura
    figure(1);
    plot(datos(p_dini:p_dfin),d(i,:),'LineWidth',2);
    hold('on');
    leyenda{i} = sprintf('DP=%.1f%%P',dp(i)*100.0);
end
hold('off');
legend(leyenda);
grid('on');
xlabel('Días usados en el ajuste desde el inicio de la erupción');
ylabel('Duración estimada de la erupción');
if pesos==0
    titulo = sprintf('Ajuste ponderado: no\n');
else
    titulo = sprintf('Ajuste ponderado: sí\n');
end
if length(saltos)==0
    titulo = sprintf('%sSaltos: no\n',titulo);
else
    titulo = sprintf(['%sSaltos: t=[',repmat(' %d,',1,length(saltos)),']'],...
                     titulo,saltos);
end
title(titulo);
%Ficheros de resultados
idf = fopen('duraciones_estimadas.txt','wb');
fprintf(idf,['%%Días ',repmat('  DP=%2d%%',1,length(dp)),...
             repmat(' Sigma_DP=%2d%%',1,length(dp)),'\n'],dp*100.0,dp*100.0);
fprintf(idf,['%3d   ',repmat('%8.3f',1,length(dp)),...
             repmat('%13.3f',1,length(dp)),'\n'],...
        [datos(p_dini:p_dfin,1) d' dt']');
fclose(idf);
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
