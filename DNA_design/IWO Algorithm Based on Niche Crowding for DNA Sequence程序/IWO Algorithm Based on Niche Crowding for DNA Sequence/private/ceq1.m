function IntValue = ceq1?( IntDNAx,IntDNAy,i)
%�ӵ�i��λ�ÿ�ʼ����Ϊc��IntDNAx��IntDNAy�Ĺ����Ӵ��ĳ���;IN:IntDNAx��IntDNAy,��������DNA,i,��ʼλ��;OUT:IntValue,������ִ�����
common=size(find(abs(IntDNAx-IntDNAy)==0),2);
abs(IntDNAx-IntDNAy)==0
common=find(abs(IntDNAx-IntDNAy)==0);
IntValue=common;
end


