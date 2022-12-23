function IntValue= s_con( IntDNAx,IntDNAy,S_con)
%Similarity函数第二项，S_con为用户指定常数;S_con为1..Length(DNA)
global S_CON;
S_CON = 6;
if nargin==2
    S_con=S_CON;
end
IntValue=0;
l=size(IntDNAx,2);
for i=1:l
    IntValue=IntValue+T(ceq(IntDNAx,IntDNAy,i),S_con);
end
end

