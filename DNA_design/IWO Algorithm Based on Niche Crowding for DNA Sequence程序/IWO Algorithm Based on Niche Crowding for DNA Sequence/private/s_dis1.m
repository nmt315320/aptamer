function IntValue = s_dis1(IntDNAx,IntDNAy,S_dis)
%Similarity������һ��(�Ż�)��S_disΪ�û�ָ������;S_disΪ0..1
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

