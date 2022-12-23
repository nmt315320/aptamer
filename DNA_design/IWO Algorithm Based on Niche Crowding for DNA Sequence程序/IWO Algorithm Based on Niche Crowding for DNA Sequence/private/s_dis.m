function IntValue = s_dis(IntDNAx,IntDNAy,S_dis)
%Similarity函数第一项，S_dis为用户指定常数;S_dis为0..1
global S_DIS;
S_DIS = 0.17;
if nargin==2
    S_dis=S_DIS;
end
Sigma_eq=0;
l=size(IntDNAx,2);
for i=1:l
    Sigma_eq=Sigma_eq+eqq(IntDNAx(i),IntDNAy(i));
end
temp=S_dis*length_nb(IntDNAy)/2;
IntValue=T(Sigma_eq,temp);
end

