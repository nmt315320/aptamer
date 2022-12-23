function [Hm,Sim]=HmSm(DNAs)
% 用于码字评价，计算DNA序列DNA1相对于在DNAs中序列的的Similarity与H-measure，方法与Shin2005论文相同，返回2个列向量，若DNA1与DNAs中的某个序列相同，则其Similarity为0
% DNA1的长度与DNAs的列数相同,相同序列的Similarity=0
% Similarity符合Shin2005，H-measure符合Shin2005,除c0720_Shin2005TSP_MEA外


% global Gap;
[m,l]=size(DNAs);
Gap1=round(l/4);  % 取序列的1/4为Gap，与Shin2005数据基本相同
for p=1:m
    SimValue=zeros(m,1);
    HmValue=zeros(m,1);
    DNA1=DNAs(p,:);
    ReverseDNA1=seqreverse(DNA1);%IntDNAx取反
    for j=1:m
        if sum(DNA1==DNAs(j,:))==l  % 相同的DNA序列不比较Similarity,仅比较H-measure
            for g=0:Gap1  %  取序列的1/4为Gap
                tempIntDNAy=[DNAs(j,:), zeros(1,g)+5, DNAs(j,:)];
                for i=-l+1:l-1
                    ShiftValue=shift(tempIntDNAy,i);
                    currentHm=h_dis(ReverseDNA1,ShiftValue)+h_con(ReverseDNA1,ShiftValue);
                    HmValue(j)=max(HmValue(j),currentHm);
                end
            end
        else                        % 同时比较Similarity与H-measure
            for g=0:Gap1  %  取序列的1/4为Gap
                tempIntDNAy=[DNAs(j,:), zeros(1,g)+5, DNAs(j,:)];
                for i=-l+1:l-1
                    ShiftValue=shift(tempIntDNAy,i);
                    currentSim=s_dis(DNA1,ShiftValue)+s_con(DNA1,ShiftValue);
                    currentHm=h_dis(ReverseDNA1,ShiftValue)+h_con(ReverseDNA1,ShiftValue);
                    SimValue(j)=max(SimValue(j),currentSim);
                    HmValue(j)=max(HmValue(j),currentHm);
                end
            end
        end
    end
Hm(p)=sum(HmValue);
Sim(p)=sum(SimValue);
end





