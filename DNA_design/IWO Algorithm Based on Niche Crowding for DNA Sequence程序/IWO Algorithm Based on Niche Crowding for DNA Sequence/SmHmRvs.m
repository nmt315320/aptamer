function [SimValue, HmValue] = SmHmRvs( IntDNAx,IntDNAy )
%����IntDNAx��IntDNAySimilarity��H-���,�ո�'-'��5����,�Ա����Գ�,�����������㣬�޸���Shin2005�о��庯������
%Similarity����Shin2005��H-measure����Shin2005,��c0720_Shin2005TSP_MEA��
%global S_DIS S_CON  H_DIS H_CON
S_DIS=0.085;S_CON=6; H_DIS=0.085; H_CON=6;
SimValue=0;
HmValue=0;
l=size(IntDNAx,2);
RIntDNAx=flip(IntDNAx);
MatchExpSim=strcat('1{',num2str(S_CON+1),',}'); % ƥ���������ȴ���S_CON��1�ַ�
MatchExpHm=strcat('1{',num2str(H_CON+1),',}'); % ƥ���������ȴ���S_CON��1�ַ�
for g=0:5
    GapIntDNAy=[IntDNAy, zeros(1,g)+5, IntDNAy];
    for i=1:l+g+1
        EqBaseVec=GapIntDNAy(i:i+l-1)==IntDNAx;% ��ͬ���01����
        EqBaseStr=strrep(num2str(EqBaseVec),' ','');
        CurrentS_dis=sum(EqBaseVec);% ��ǰS_dis
        if CurrentS_dis>S_DIS*l
            S_dis=CurrentS_dis;
        else
            S_dis=0;
        end
        [SmMatchStar,SmMatchEnd]=regexp(EqBaseStr,MatchExpSim);
        S_con=sum(SmMatchEnd-SmMatchStar+1);
        currentSim=S_dis+S_con;
        SimValue=max(currentSim,SimValue);

        BpBaseVec=((GapIntDNAy(i:i+l-1)+RIntDNAx)==5);  % �������01����
        BpBaseStr=strrep(num2str(BpBaseVec),' ','');
        CurrentH_dis=sum(BpBaseVec);% ��ǰH_dis
        if CurrentH_dis>H_DIS*l
            H_dis=CurrentH_dis;
        else
            H_dis=0;
        end
        
       % MatchExp     
        [HmMatchStar,HmMatchEnd]=regexp(BpBaseStr,MatchExpHm);
       % 'H_con'
        H_con=sum(HmMatchEnd-HmMatchStar+1);
        currentHm=H_dis+H_con;
        HmValue=max(currentHm,HmValue);
    end
end
end
