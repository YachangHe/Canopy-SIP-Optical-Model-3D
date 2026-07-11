function  bg =get_HSF_go( par,SZA,SAA,VZA,VAA,Ps_dir_go,Pv_dir_go,z)
%This routine calculates the bidirectional between-crown gap probability 
%KG = Pw(omega0)Pw(omega) + Y*sqrt(Pw(omega0)Pw(omega)(1-Pw(omega0))(1-Pw(omega0)))
%input:
%par:crown size
%geometry in degree!
%z:crown center height  (Height - lmax + 1/2*Crowndeepth)
%output:
%bg:bidirectional between-crown gap probability KG

%%
SZA=deg2rad(SZA);SAA=deg2rad(SAA);
VZA=deg2rad(VZA);VAA=deg2rad(VAA);
mu0=cos(SZA);
muv=cos(VZA);
f1=sqrt(Ps_dir_go*Pv_dir_go*(1-Ps_dir_go)*(1-Pv_dir_go));     % 协方差；对应文章(7) 根号下的式子
%calculate delta:
phi = VAA-SAA;
if(phi < 0 && phi > - 360)
   phi = phi + 360;
elseif(phi > 360 && phi < 720)
   phi = phi - 360;
else
   phi = phi;
end
cosgamma=cos(SZA)*cos(VZA)+sin(SZA)*sin(VZA)*cos(VAA-SAA);
delta=sqrt(1/(cos(SZA)^2)+1/(cos(VZA)^2)-2*cosgamma/(cos(SZA)*cos(VZA)));     % 公式(9): angular distance measure
if delta<0.00001
    delta=0.00001;
end
Y=exp(-delta/par*z);
bg  = f1 * Y;

end

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     