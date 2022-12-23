function IntValue = s_dis1(IntDNAx,IntDNAy,S_dis)
%Similarity函数第一项(优化)，S_dis为用户指定常数;S_dis为0..1
global S_DIS;
if nargin==2
    S_dis=S_DIS;
end
IntValue=0;
Sigma_eq=size(find(abs(IntDNAx-IntDNAy)==0),2);
temp=S_dis*length_nb(IntDNAy)/2;
if Sigma_eq>temp
    IntValue=Sigma_eq;
end
end

