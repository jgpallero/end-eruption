%AJUSTE DE S0, SALTOS Y TAU
%Carga de datos
%Fichero de tres columnas: número de día, señal a ajustar, desviación típica
datos = load('acortamiento-mazo-lp01.txt');
% datos = load('volumenes-material.txt');
%Días (incluidos) entre los que se usarán observaciones para el ajuste
%Se usarán>=diaIni y <=diaFin
diaIni = 1;
diaFin = 70;
%Saltos (números de los días, todos los días >= están en el salto)
% saltos = [0 45 49 62];
saltos = [62];
% saltos = [];
%Uso de ponderación: 0/1 -> no/sí
pesos = 1;
%Criterios de parada (norm(dx)/norm(x) y número máximo de iteraciones)
parada = [0.0001 20];
%Detección de errores groseros (eg==0 o eg>=1 no los busca, valor en (0,1)
%quiere decir el nivel de significación del test de Pope y un número negativo es
%el residuo máximo)
% eg = 0.01;
eg = -5.0;
% eg = 0;
%Días (incluidos) desde y hasta que se usarán datos para dibujar
ddi = 0;
ddf = 160;
%Intervalo de confianza para dibujar los márgenes de error del modelo
ic = 0.99;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Posiciones de los días a usar
p_dias_ini = find(datos(:,1)>=diaIni);
if length(p_dias_ini)==0
    error(['El día indicado en la variable ''diaIni'' es mayor que todos ',...
           'los días almacenados en el fichero']);
end
p_dias_ini = p_dias_ini(1);
p_dias_fin = find(datos(:,1)<=diaFin);
if length(p_dias_fin)==0
    error(['El día indicado en la variable ''diaFin'' es menor que todos ',...
           'los días almacenados en el fichero']);
end
p_dias_fin = p_dias_fin(end);
%Posiciones de los días para dibujar
p_dd_ini = find(datos(:,1)>=ddi);
if length(p_dd_ini)==0
    error(['El día indicado en la variable ''ddi'' es mayor que todos los ',...
           'días almacenados en el fichero']);
end
p_dd_ini = p_dd_ini(1);
p_dd_fin = find(datos(:,1)<=ddf);
if length(p_dd_fin)==0
    error(['El día indicado en la variable ''ddf'' es menor que todos los ',...
           'días almacenados en el fichero']);
end
p_dd_fin = p_dd_fin(end);
%Datos de trabajo para el ajuste y para el dibujo
datos_a = datos(p_dias_ini:p_dias_fin,:);
datos_d = datos(p_dd_ini:p_dd_fin,:);
%Realizo el ajuste
[s0,tau,c,Exx,iter,du] = AjustaModExp(datos_a,saltos,pesos,eg,parada);
if iter==parada(2)
    error('Se ha llegado al número máximo de iteraciones para el ajuste');
end
%Desviaciones típicas de los parámetros
sx = sqrt(diag(Exx));
%Calculamos con el modelo
sc = s0*(1.0-exp(-datos_d(:,1)./tau));
for i=1:size(c,2)
    pos = datos_d(:,1)>=c(1,i);
    sc(pos) = sc(pos)+c(2,i);
end
%Propagación de varianzas y covarianzas
sxx = zeros(size(sc));
J = zeros(1,2+size(c,2));
datos_dd = [];
for i=1:length(sc)
    %Matriz jacobiana
    J(1) = 1.0-exp(-datos_d(i,1)./tau);
    J(2) = -s0*datos_d(i,1)./(tau^2).*exp(-datos_d(i,1)./tau);
    for j=1:size(c,2)
        if datos_d(i,1)>=c(1,j)
            J(2+j) = 1.0;
        else
            J(2+j) = 0.0;
        end
    end
    sxx(i) = sqrt(J*Exx*J');
    %Compruebo si el dato no está entre los usados en el cálculo
    if ~sum(datos_d(i,1)==datos_a(:,1))
        datos_dd = [datos_dd;datos_d(i,:)];
    end
end
%Figura
errorbar(datos_a(du,1),datos_a(du,2),datos_a(du,3),'ob');
hold('on');
if sum(~du)
    errorbar(datos_a(~du,1),datos_a(~du,2),datos_a(~du,3),'^m');
end
if size(datos_dd,1)>0
    errorbar(datos_dd(:,1),datos_dd(:,2),datos_dd(:,3),'og');
end
se = abs(norminv((1.0-ic)/2))*sxx;
plot(datos_d(:,1),sc,'r','LineWidth',2);
plot(datos_d(:,1),sc+se,'r--','LineWidth',2);
plot(datos_d(:,1),sc-se,'r--','LineWidth',2);
hold('off');
xlim([min(datos_d(:,1)) max(datos_d(:,1))]);
xlabel('Días desde el inicio de la erupción');
ylabel('Señal observada');
titulo = sprintf('Ajuste ponderado: ');
if pesos==0
    titulo = sprintf('%sno, ',titulo);
else
    titulo = sprintf('%ssí, ',titulo);
end
titulo = sprintf('%sdetección errores groseros: ',titulo);
if (eg==0.0)||(eg>=1.0)
    titulo = sprintf('%sno\n',titulo);
elseif eg<0.0
    titulo = sprintf('%sresiduo>%.3f\n',titulo,abs(eg));
else
    titulo = sprintf('%sPope (alfa=%.3f)\n',titulo,abs(eg));
end
titulo = sprintf('%sDías usados en el ajuste: %d',titulo,length(du));
if sum(~du)
    titulo = sprintf('%s, días eliminados: %d',titulo,sum(~du));
    titulo = sprintf('%s, días finalmente utilizados: %d\n',titulo,sum(du));
else
    titulo = sprintf('%s\n',titulo);
end
if size(c,2)==0
    titulo = sprintf('%sParámetros c(t-t_k): no\n',titulo);
else
    titulo = sprintf(['%sParámetros c(t-t_k): t_k=[',...
                      repmat(' %d,',1,size(c,2)),'] --> c(t-t_k)=['],...
                     titulo,c(1,:));
    titulo = sprintf(['%s',repmat(' %.4f,',1,size(c,2)),'] +/- ['],...
                      titulo,c(2,:));
    titulo = sprintf(['%s',repmat(' %.4f,',1,size(c,2)),'] (sigma)\n'],...
                     titulo,sx(3:end));
end
titulo = sprintf('%ss0=%.4f +/- %.4f (sigma), tau=%.4f +/- %.4f (sigma)',...
                 titulo,s0,sx(1),tau,sx(2));
title(titulo);
orden_leyenda = 'legend(''Datos usados en el ajuste''';
if sum(~du)
    orden_leyenda = [orden_leyenda,',''Datos rechazados'''];
end
if size(datos_dd,1)>0
    orden_leyenda = [orden_leyenda,',''Datos no usados'''];
end
orden_leyenda = [orden_leyenda,');'];
eval(orden_leyenda);
grid('on');
print('modelo_ajustado.png','-dpng');
%Ficheros de resultados
idf = fopen('puntos_usados.txt','wb');
fprintf(idf,'%%Observaciones originales\n');
fprintf(idf,'%%Puntos usados en el ajuste (errores groseros excluidos)\n');
fprintf(idf,'%%Día  Acort.  (mm)   Sigma (mm)\n');
fprintf(idf,['%3d ',repmat(' %14.7E',1,size(datos_a,2)-1),'\n'],datos_a(du,:)');
fclose(idf);
idf = fopen('puntos_modelo.txt','wb');
fprintf(idf,'%%Puntos calculados con el modelo\n');
fprintf(idf,'%%El intervalo de confianza es %.3f\n',ic);
fprintf(idf,'%%Día  Acort.  (mm)   Sigma (mm)\n');
fprintf(idf,['%3d  %14.7E %14.7E\n'],[datos_d(:,1) sc se]');
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
