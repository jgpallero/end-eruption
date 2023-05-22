%Load a text file (columns separated by spaces):
%Number_of_day  Baseline_shortening  Stdv_Baseline_Shortening
datos = load('acortamiento-mazo-lp01.txt');
%Día (incluido) hasta el que se usarán observaciones para el ajuste)
dias = 70;
%Saltos (número del día), uso de ponderación (0/1 -> no/sí) y criterio de parada
saltos = [62];
% saltos = [];
pesos = 1;
parada = [0.0001 10];
%Detección de errores groseros (eg==0 o eg>=1 no los busca, valor en (0,1)
%quiere decir el nivel de significación del test de Pope y un número negativo es
%el residuo máximo)
% eg = 0.01;
% eg = -5.0;
eg = 0;
%Día (incluido) hasta el que se usarán datos para dibujar
dd = 125;
%Intervalo de confianza para dibujar los márgenes de error del modelo
ic = 0.99;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Posiciones de los días a usar y oara dibujar
p_dias = find(datos(:,1)==dias);
if length(p_dias)==0
    error('El día indicado en la variable ''dias'' no existe en el fichero');
end
p_dd = find(datos(:,1)==dd);
if length(p_dd)==0
    error('El día indicado en la variable ''dd'' no existe en el fichero');
end
if p_dd<p_dias
    error(['El valor de la variable ''dias'' ha de ser mayor que el de la ',...
           'variable''dd''']);
end
%Datos de trabajo
datos = datos(1:p_dd,:);
%Realizo el ajuste
[s0,tau,c,Exx,iter,du] = AjustaModExp(datos(1:p_dias,:),saltos,pesos,eg,parada);
if iter==parada(2)
    error('Se ha llegado al número máximo de iteraciones para el ajuste');
end
%Vector de datos usados y no usados
ndu = sum(du);
dumax = max(datos(du,1));
du = [du;logical(ones(size(datos,1)-length(du),1))];
%Desviaciones típicas de los parámetros
sx = sqrt(diag(Exx));
%Calculamos con el modelo
sc = s0*(1.0-exp(-datos(:,1)./tau));
for i=1:size(c,2)
    pos = datos(:,1)>=c(1,i);
    sc(pos) = sc(pos)+c(2,i);
end
%Propagación de varianzas y covarianzas
sxx = zeros(size(sc));
J = zeros(1,2+size(c,2));
for i=1:length(sc)
    %Matriz jacobiana
    J(1) = 1.0-exp(-datos(i,1)./tau);
    J(2) = -s0*datos(i,1)./(tau^2).*exp(-datos(i,1)./tau);
    for j=1:size(c,2)
        if datos(i,1)>=c(1,j)
            J(2+j) = 1.0;
        else
            J(2+j) = 0.0;
        end
    end
    sxx(i) = sqrt(J*Exx*J');
end
%Figura
errorbar(datos(du,1),datos(du,2),datos(du,3),'ob');
hold('on');
if sum(~du)
    errorbar(datos(~du,1),datos(~du,2),datos(~du,3),'^m');
end
se = abs(norminv((1.0-ic)/2))*sxx;
plot(datos(:,1),sc,'r','LineWidth',2);
plot(datos(:,1),sc+se,'r--','LineWidth',2);
plot(datos(:,1),sc-se,'r--','LineWidth',2);
hold('off');
xlim([min(datos(:,1)) max(datos(:,1))]);
xlabel('Días desde el inicio de la erupcion');
ylabel('Acortamiento de la línea');
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
titulo = sprintf('%sDías usados en el ajuste: %d',titulo,p_dias);
if sum(~du)
    titulo = sprintf('%s, días eliminados: %d',titulo,sum(~du));
    titulo = sprintf('%s, días finalmente utilizados: %d\n',titulo,ndu);
else
    titulo = sprintf('%s\n',titulo);
end
if size(c,2)==0
    titulo = sprintf('%sParámetros c(t-t_k): no\n',titulo);
else
    titulo = sprintf(['%sParámetros c(t-t_k): t_k=[',...
                      repmat(' %d,',1,size(c,2)),'] --> c(t-t_k)=['],...
                     titulo,c(1,:));
    titulo = sprintf(['%s',repmat(' %.3f,',1,size(c,2)),'] +/- ['],...
                      titulo,c(2,:));
    titulo = sprintf(['%s',repmat(' %.3f,',1,size(c,2)),'] (sigma)\n'],...
                     titulo,sx(3:end));
end
titulo = sprintf('%ss0=%.3f +/- %.3f (sigma), tau=%.3f +/- %.3f (sigma)',...
                 titulo,s0,sx(1),tau,sx(2));
title(titulo);
grid('on');
%Ficheros de resultados
idf = fopen('puntos_usados.txt','wb');
fprintf(idf,'%%Observaciones originales\n');
fprintf(idf,'%%En el ajuste se han usado los puntos hasta el día %d\n',dumax);
fprintf(idf,'%%Día  Acort.  (mm)   Sigma (mm)\n');
fprintf(idf,['%3d ',repmat(' %14.7E',1,size(datos,2)-1),'\n'],datos(du,:)');
fclose(idf);
idf = fopen('puntos_modelo.txt','wb');
fprintf(idf,'%%Puntos calculados con el modelo hasta el día %d\n',dd);
fprintf(idf,'%%El intervalo de confianza es %.3f\n',ic);
fprintf(idf,'%%Día  Acort.  (mm)   Sigma (mm)\n');
fprintf(idf,['%3d  %14.7E %14.7E\n'],[datos(:,1) sc se]');
fclose(idf);
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
