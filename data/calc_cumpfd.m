load("fy_U5.txt");

cumpdf=zeros(length(fy_U5),1);
x=[1:119];
integ=0;

for i=1:length(cumpdf)

integ=integ+fy_U5(i,2);
cumpdf(i)=integ;

end

plot(cumpdf./200,x);

fich=fopen("proba_fy_U5_pourintegration_dansopendec.txt","w");


fprintf(fich,"%e, ",cumpdf./200);

fclose(fich);

