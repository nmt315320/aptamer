function [Hm,Sim]=HmSm(DNAs)
% �����������ۣ�����DNA����DNA1�������DNAs�����еĵ�Similarity��H-measure��������Shin2005������ͬ������2������������DNA1��DNAs�е�ĳ��������ͬ������SimilarityΪ0
% DNA1�ĳ�����DNAs��������ͬ,��ͬ���е�Similarity=0
% Similarity����Shin2005��H-measure����Shin2005,��c0720_Shin2005TSP_MEA��


% global Gap;
[m,l]=size(DNAs);
Gap1=round(l/4);  % ȡ���е�1/4ΪGap����Shin2005���ݻ�����ͬ
for p=1:m
    SimValue=zeros(m,1);
    HmValue=zeros(m,1);
    DNA1=DNAs(p,:);
    ReverseDNA1=seqreverse(DNA1);%IntDNAxȡ��
    for j=1:m
        if sum(DNA1==DNAs(j,:))==l  % ��ͬ��DNA���в��Ƚ�Similarity,���Ƚ�H-measure
            for g=0:Gap1  %  ȡ���е�1/4ΪGap
                tempIntDNAy=[DNAs(j,:), zeros(1,g)+5, DNAs(j,:)];
                for i=-l+1:l-1
                    ShiftValue=shift(tempIntDNAy,i);
                    currentHm=h_dis(ReverseDNA1,ShiftValue)+h_con(ReverseDNA1,ShiftValue);
                    HmValue(j)=max(HmValue(j),currentHm);
                end
            end
        else                        % ͬʱ�Ƚ�Similarity��H-measure
            for g=0:Gap1  %  ȡ���е�1/4ΪGap
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





