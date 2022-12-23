function IntValue= h_dis( IntDNAx,IntDNAy,H_dis)
%H_measure函数第一项，H_dis为用户指定常数;H_dis为0..1
global H_DIS;
H_DIS = 0.17;
if nargin==2
    H_dis=H_DIS;
end
Sigma_bp=0;
l=size(IntDNAx,2);
for i=1:l
    Sigma_bp=Sigma_bp+bp(IntDNAx(i),IntDNAy(i));
end
temp=H_dis*length_nb(IntDNAy)/2;
IntValue=T(Sigma_bp,temp);
end

