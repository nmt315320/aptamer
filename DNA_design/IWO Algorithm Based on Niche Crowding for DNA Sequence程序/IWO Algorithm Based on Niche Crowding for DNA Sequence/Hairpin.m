function Values = Hairpin( DNAs )
% 计算DNAs(整数编码矩阵)的发卡数，返回Values列向量
% function IntValue = Hairpin( IntDNA, R_min,P_min )
%依Shin2005,计算IntDNA的发卡数
global R_MIN P_MIN  PINLEN
R_MIN=6;
P_MIN=6;
PINLEN=3;
if nargin==1
    R_min=R_MIN;
    P_min=P_MIN;
end
[m,l]=size(DNAs);
Values=zeros(m,1);
for k=1:m
    IntDNA=DNAs(k,:);
    for p=P_min:(l-R_min)/2
        for r=R_min:(l-2*p)
            for i=0:(l-2*p-r)
                Sigma_bp=0;
                for j=0:pinlen(p,r,i,IntDNA)-1%%j=1..pinlen-1   %%j=1...pinlen
                    Sigma_bp=Sigma_bp+bp(IntDNA(p+i-j),IntDNA(p+i+r+1+j));      %%Sigma_bp=Sigma_bp+bp(IntDNA(p+i-j),IntDNA(p+i+r+j));
                end
                Values(k)=Values(k)+T(Sigma_bp,pinlen(p,r,i,IntDNA)/ PINLEN);
            end
        end
    end
end
end

function minValue = pinlen( p,r,i,IntDNA )
%辅助函数，返回min(p+i,l-(p+i+r))
l=size(IntDNA,2);
minValue=min(p+i,l-p-i-r);
end
