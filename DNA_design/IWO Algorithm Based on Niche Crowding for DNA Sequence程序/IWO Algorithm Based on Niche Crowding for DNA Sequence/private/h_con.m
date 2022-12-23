function IntValue = h_con( IntDNAx,IntDNAy,H_con )
%H_measure函数第二项，H_con为用户指定常数;H_con为1..Length(DNA)
global H_CON;
H_CON=6;
if nargin==2
    H_con=H_CON;
end
IntValue=0;
l=size(IntDNAx,2);
for i=1:l
    IntValue=IntValue+T(cbp(IntDNAx,IntDNAy,i),H_con);
end
end

