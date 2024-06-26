%binarychanneldetection.m
function [datatx,datarx,datadetected,...
    bayesruleselected,...
    pl,ph,val_pl,val_ph,...
    pd,gamma,pd_obtained_beforedetection,...
    pfa_obtained_beforedetection,...
    pd_obtained_afterdetection,...
    pfa_obtained_afterdetection]...
   = binarychanneldetection(x,y,z,len,p0,pfa)
datatx=[];
temp=round(rand(1,len));
for i=1:1:len
r=rand;
u=[p0 1];
test=u-r;
if(test(1)>0)
datatx=[datatx 0];
else
datatx=[datatx 1];
end
end
%x->vector [p(0/0) p(1/1)]
%y=[c01 c10];
%z->1:Bayes technique.
%z->2:Mini-max technique.
%z->3:Neyman-pearson technique needs probabily of false alarm pfa.
%p0->prior probability
%datatx->data transmitted
%datarx->data received
%bayesruleselected
%pl-lower level probability used in mini-max technique
%ph-upper level probability used in mini-max technique
%val_pl-bayes cost at pl
%val_ph-bayes cost at ph
%pd-probability of detection computed for four choices of Neyman-pearson
%technique.
%gamma-probabilities used in four choices of Neyman-pearson
%technique.
%pd_obtained-probability of detection actually obtained in the simulation.
%pfa_obtained-probability of false alarm actually obtained in the
%simulation.
p00=x(1);
p10=1-x(1);
p11=x(2);
p01=1-x(2);
datarx=[];
for i=1:1:length(datatx)
if(datatx(i)==0)
r=rand;
u=[p00 1];
test=u-r;
if(test(1)>0)
datarx=[datarx 0];
else
datarx=[datarx 1];
end
else
r=rand;
u=[p01 1];
test=u-r;
if(test(1)>0)
datarx=[datarx 0];
else
datarx=[datarx 1];
end
end
end
switch(z)
case 1
pl=[];ph=[];val_pl=[];val_ph=[];pd=[];gamma=[];
r1cost=y(2)*p0*p10+y(1)*(1-p0)*p01;
r2cost=y(2)*p0*p00+y(1)*(1-p0)*p11;
r3cost=1-p0;
r4cost=p0;110 3 Detection Theory and Estimation Theory for Wireless Communication
[u,v]=min([r1cost r2cost r3cost r4cost]);
bayesruleselected=v;
switch(v)
case 1
datadetected=datarx;
case 2
datadetected=[];
for i=1:1:length(datarx)
if(datarx(i)==0)
datadetected=[datadetected 1];
else
datadetected=[datadetected 0];
end
end
case 3
datadetected=zeros(1,length(datarx));
case 4
datadetected=ones(1,length(datarx));
end
case 2
%Mini-max technique
pd=[];gamma=[];
bayesruleselected=[];
pl=min((p11/(1-p00+p11)),(p01/(1-p10+p01)));
ph=max((p11/(1-p00+p11)),(p01/(1-p10+p01)));
r1costpl=y(2)*pl*p10+y(1)*(1-pl)*p01;
r2costpl=y(2)*pl*p00+y(1)*(1-pl)*p11;
r3costpl=1-pl;
r4costpl=pl;
[val_pl bayes_pl]=min([r1costpl r2costpl r3costpl r4costpl]);
r1costph=y(2)*ph*p10+y(1)*(1-ph)*p01;
r2costph=y(2)*ph*p00+y(1)*(1-ph)*p11;
r3costph=1-ph;
r4costph=ph;
[val_ph bayes_ph]=min([r1costph r2costph r3costph r4costph]);
[minimaxval minimaxpos]=max([val_pl val_ph]);
%Bayes rule corresponding to the p0 ranging between pl and ph
pchoose=(pl+ph)/2;
r1costpchoose=y(2)*pchoose*p10+y(1)*(1-pchoose)*p01;
r2costpchoose=y(2)*pchoose*p00+y(1)*(1-pchoose)*p11;
r3costpchoose=1-pchoose;
r4costpchoose=pchoose;
[val_pchoose bayes_pchoose]=min([r1costpchoose r2costpchoose ...
r3costpchoose r4costpchoose]);
if((ph==pl))
minimaxpos=3;
end
switch(minimaxpos)
case 1
%Randomized decision rule between the region 1 and 2.
%Randomized decision rule between the rule 4 with
%probability rho and rule described
%by the variable bayes_pchoose with probability 1-rho
temp=(val_pl-val_ph+eps)/(ph-pl+eps);
rho=temp/(temp+1)
datadetected=[];
for w=1:1:len
t1=[rho 1];
u=t1-rand;
if(u(1)>0)3.1 Detection Theory for Binary Signal Transmission 111
datadetected=[datadetected 1];
else
switch(bayes_pchoose)
case 1
datadetected=[datadetected datarx(w)];
case 2
if(datarx(w)==0)
datadetected=[datadetected 1];
else
datadetected=[datadetected 0];
end
end
end
end
case 2
%Randomized decision rule between the region 2 and 3.
%Randomized decision rule between the rule described
%by the variable bayes_pchoose and rule 3
temp=(val_ph-val_pl+eps)/(ph-pl+eps);
rho=1/(1+temp);
datadetected=[];
for w=1:1:len
t1=[rho 1];
u=t1-rand;
if(u(1)>0)
switch(bayes_pchoose)
case 1
datadetected=[datadetected datarx(w)];
case 2
if(datarx(w)==0)
datadetected=[datadetected 1];
else
datadetected=[datadetected 0];
end
end
else
datadetected=[datadetected 0];
end
end
case 3
datadetected=[];
for w=1:1:len
t1=[0.5 1];
u=t1-rand;
if(u(1)>0)
datadetected=[datadetected 1];
else
datadetected=[datadetected 0];
end
end
end
case 3
%Neyman-pearson technique
pl=[];ph=[];val_pl=[];val_ph=[];pd=[];
%when the received signal is 1, decide it as in favour of 1 with
%probability 1. when the received signal is 0, decide it as in favour of 1112 3 Detection Theory and Estimation Theory for Wireless Communication
%with probability gamma(1).
gamma(1)=(pfa-p10)/p00;
pd(1)=p11+gamma(1)*p01;
%when the received signal is 1, decide it as in favour of 1 with
%probability gamma(2). when the received signal is 0, decide it
%as in favour of 1 with probability 0.
gamma(2)=pfa/p10;
pd(2)=gamma(2)*p11;
%when the received signal is 0, decide it as in favour of 1 with
%probability 1. when the received signal is 1, decide it as in favour of 1
%with probability gamma(3).
gamma(3)=(pfa-p00)/p10;
pd(3)=p00+gamma(3)*p10;
%when the received signal is 0, decide it as in favour of 1 with
%probability gamma(4). when the received signal is 1, decide
%it as in favour of 1 with probability 0.
gamma(4)=pfa/p00;
pd(4)=gamma(4)*p00+0*p10;
validpd=[];
for w=1:1:4
if(gamma(w)>=0)
validpd=[validpd pd(w)];
else
validpd=[validpd -1];
end
end
[p,q]=max(validpd);
switch(q)
case 1
datadetected=[];
for c=1:1:length(datarx)
if(datarx(c)==1)
datadetected=[datadetected 1];
else
t1=[gamma(1) 1];
u=t1-rand;
if(u(1)>0)
datadetected=[datadetected 1];
else
datadetected=[datadetected 0];
end
end
end
case 2
datadetected=[];
for c=1:1:length(datarx)
if(datarx(c)==0)
datadetected=[datadetected 0];
else
t1=[gamma(2) 1];
u=t1-rand;
if(u(1)>0)
datadetected=[datadetected 1];
else
datadetected=[datadetected 0];
end
end
end
case 33.1 Detection Theory for Binary Signal Transmission 113
datadetected=[];
for c=1:1:length(datarx)
if(datarx(c)==0)
datadetected=[datadetected 1];
else
t1=[gamma(3) 1];
u=t1-rand;
if(u(1)>0)
datadetected=[datadetected 1];
else
datadetected=[datadetected 0];
end
end
end
case 4
datadetected=[];
for c=1:1:length(datarx)
if(datarx(c)==1)
datadetected=[datadetected 0];
else
t1=[gamma(4) 1];
u=t1-rand;
if(u(1)>0)
datadetected=[datadetected 1];
else
datadetected=[datadetected 0];
end
end
end
end
bayesruleselected=[];
pl=[];ph=[];bayes_ph=[];bayes_pl=[];
end
pd_obtained_afterdetection=(length(find((datadetected-datatx)==0))...
/length(datatx));
pd_obtained_beforedetection=(length(find((datarx-datatx)==0))...
/length(datatx));
pfa_obtained_afterdetection=(length(find((datadetected-datatx)==1))...
/length(datatx));
pfa_obtained_beforedetection=(length(find((datarx-datatx)==1))...
/length(datatx));
figure(1)
subplot(3,1,1)
plot(datatx,’r’)
subplot(3,1,2)
plot(datarx,’g’)
subplot(3,1,3)
plot(datadetected,’b’)
bayes.m
close all
y=[1 1];
z=1;
len=1000;
p0=0.5;
pfa=[];
CLUSTER1=[];
CLUSTER2=[];
CLUSTER3=[];
CLUSTER4=[];114 3 Detection Theory and Estimation Theory for Wireless Communication
for j=1:1:1000
x=[rand rand];
[datatx,datarx,datadetected,bayesruleselected,pl,ph,val_pl,val_ph,...
pd,gamma,pd_obtained,pfa_obtained]...
=binarychanneldetection(x,y,z,len,p0,pfa)
v=bayesruleselected;
if(v==1)
CLUSTER1=[CLUSTER1;x];
elseif(v==2)
CLUSTER2=[CLUSTER2;x];
elseif(v==3)
CLUSTER3=[CLUSTER3;x];
else
CLUSTER4=[CLUSTER4;x];
end
end
p0=0.9;
[CLUSTER1, CLUSTER2, CLUSTER3, CLUSTER4]= bayes(p0);
figure(10)
if(isempty(CLUSTER1)==0)
plot(CLUSTER1(:,1),CLUSTER1(:,2),’r*’)
end
hold on
if(isempty(CLUSTER2)==0)
plot(CLUSTER2(:,1),CLUSTER2(:,2),’go’)
end
if(isempty(CLUSTER3)==0)
plot(CLUSTER3(:,1),CLUSTER3(:,2),’b+’)
end
if(isempty(CLUSTER4)==0)
plot(CLUSTER4(:,1),CLUSTER4(:,2),’k+’)
end
p0=0.7;
[CLUSTER1, CLUSTER2, CLUSTER3, CLUSTER4]= bayes(p0);
figure(11)
if(isempty(CLUSTER1)==0)
plot(CLUSTER1(:,1),CLUSTER1(:,2),’r*’)
end
hold on
if(isempty(CLUSTER2)==0)
plot(CLUSTER2(:,1),CLUSTER2(:,2),’go’)
end
if(isempty(CLUSTER3)==0)
plot(CLUSTER3(:,1),CLUSTER3(:,2),’b+’)
end
if(isempty(CLUSTER4)==0)
plot(CLUSTER4(:,1),CLUSTER4(:,2),’k+’)
end
p0=0.5;
[CLUSTER1, CLUSTER2, CLUSTER3, CLUSTER4]= bayes(p0);
[CLUSTER1, CLUSTER2, CLUSTER3, CLUSTER4]= bayes(p0);
figure(12)
if(isempty(CLUSTER1)==0)
plot(CLUSTER1(:,1),CLUSTER1(:,2),’r*’)
end
hold on
if(isempty(CLUSTER2)==0)
plot(CLUSTER2(:,1),CLUSTER2(:,2),’go’)
end3.1 Detection Theory for Binary Signal Transmission 115
if(isempty(CLUSTER3)==0)
plot(CLUSTER3(:,1),CLUSTER3(:,2),’b+’)
end
if(isempty(CLUSTER4)==0)
plot(CLUSTER4(:,1),CLUSTER4(:,2),’k+’)
end
p0=0.3;
[CLUSTER1, CLUSTER2, CLUSTER3, CLUSTER4]= bayes(p0);
figure(13)
if(isempty(CLUSTER1)==0)
plot(CLUSTER1(:,1),CLUSTER1(:,2),’r*’)
end
hold on
if(isempty(CLUSTER2)==0)
plot(CLUSTER2(:,1),CLUSTER2(:,2),’go’)
end
if(isempty(CLUSTER3)==0)
plot(CLUSTER3(:,1),CLUSTER3(:,2),’b+’)
end
if(isempty(CLUSTER4)==0)
plot(CLUSTER4(:,1),CLUSTER4(:,2),’k+’)
end
p0=0.1;
[CLUSTER1, CLUSTER2, CLUSTER3, CLUSTER4]= bayes(p0);
figure(14)
if(isempty(CLUSTER1)==0)
plot(CLUSTER1(:,1),CLUSTER1(:,2),’r*’)
end
hold on
if(isempty(CLUSTER2)==0)
plot(CLUSTER2(:,1),CLUSTER2(:,2),’go’)
end
if(isempty(CLUSTER3)==0)
plot(CLUSTER3(:,1),CLUSTER3(:,2),’b+’)
end
if(isempty(CLUSTER4)==0)
plot(CLUSTER4(:,1),CLUSTER4(:,2),’k+’)
end
%minimax.m
close all
y=[1 1];
z=2;
len=1000;
priorprob{1}=[0.8 0.7];
priorprob{2}=[0.4 0.8];
priorprob{3}=[0.7 0.2];
priorprob{4}=[0.4 0.2];
pfa=[];
for i=1:1:4
x=priorprob{i};
pdcol_before=[];
pfacol_before=[];
pdcol_after=[];
pfacol_after=[];
for p0=[0.1:1/1000:0.9]
[datatx,datarx,datadetected,bayesruleselected,pl,ph,val_pl,val_ph,...
pd,gamma,pd_obtained_beforedetection,pfa_obtained_beforedetection,...116 3 Detection Theory and Estimation Theory for Wireless Communication
pd_obtained_afterdetection,pfa_obtained_afterdetection]...
=binarychanneldetection(x,y,z,len,p0,pfa);
pdcol_before=[pdcol_before pd_obtained_beforedetection];
pfacol_before=[pfacol_before pfa_obtained_beforedetection];
pdcol_after=[pdcol_after pd_obtained_afterdetection];
pfacol_after=[pfacol_after pfa_obtained_afterdetection];
end
pdcol_beforefinal{i}=pdcol_before;
pfacol_beforefinal{i}=pfacol_before;
pdcol_afterfinal{i}=pdcol_after;
pfacol_afterfinal{i}=pfacol_after;
end
for i=1:1:4
figure(2*(i-1)+2)
plot(pfacol_beforefinal{i},pdcol_beforefinal{i},’*’)
figure(2*(i-1)+3)
plot(pfacol_afterfinal{i},pdcol_afterfinal{i},’*’)
end
close all
%neymanpearson.m
y=[1 1];
z=3;
pfa=0.01;
len=1000;
pdcolbeforedetection1=[];
pfacolbeforedetection1=[];
pdcolafterdetection1=[];
pfacolafterdetection1=[];
x=[0.8 0.7];
for p0=0.01:1/1000:0.99
[datatx,datarx,datadetected,bayesruleselected,pl,ph,val_pl,val_ph,...
pd,gamma,pd_obtained_beforedetection,pfa_obtained_beforedetection,...
pd_obtained_afterdetection,pfa_obtained_afterdetection]...
=binarychanneldetection(x,y,z,len,p0,pfa);
pdcolbeforedetection1=[pdcolbeforedetection1 pd_obtained_beforedetection];
pfacolbeforedetection1=[pfacolbeforedetection1 pfa_obtained_beforedetection];
pdcolafterdetection1=[pdcolafterdetection1 pd_obtained_afterdetection];
pfacolafterdetection1=[pfacolafterdetection1 pfa_obtained_afterdetection];
end
figure
plot(pfacolbeforedetection1,pdcolbeforedetection1,’*’)
figure
plot(pfacolafterdetection1,pdcolafterdetection1,’*’)
pdcolbeforedetection2=[];
pfacolbeforedetection2=[];
pdcolafterdetection2=[];
pfacolafterdetection2=[];
x=[0.4 0.8];
for p0=0.01:1/1000:0.99
[datatx,datarx,datadetected,bayesruleselected,pl,ph,val_pl,val_ph,...
pd,gamma,pd_obtained_beforedetection,pfa_obtained_beforedetection,...
pd_obtained_afterdetection,pfa_obtained_afterdetection]...
=binarychanneldetection(x,y,z,len,p0,pfa);
pdcolbeforedetection2=[pdcolbeforedetection2 pd_obtained_beforedetection];
pfacolbeforedetection2=[pfacolbeforedetection2 pfa_obtained_beforedetection];
pdcolafterdetection2=[pdcolafterdetection2 pd_obtained_afterdetection];
pfacolafterdetection2=[pfacolafterdetection2 pfa_obtained_afterdetection];
end3.1 Detection Theory for Binary Signal Transmission 117
figure
plot(pfacolbeforedetection2,pdcolbeforedetection2,’*’)
figure
plot(pfacolafterdetection2,pdcolafterdetection2,’*’)
pdcolbeforedetection3=[];
pfacolbeforedetection3=[];
pdcolafterdetection3=[];
pfacolafterdetection3=[];
x=[0.7 0.2];
for p0=0.01:1/1000:0.99
[datatx,datarx,datadetected,bayesruleselected,pl,ph,val_pl,val_ph,...
pd,gamma,pd_obtained_beforedetection,pfa_obtained_beforedetection,...
pd_obtained_afterdetection,pfa_obtained_afterdetection]...
=binarychanneldetection(x,y,z,len,p0,pfa);
pdcolbeforedetection3=[pdcolbeforedetection3 pd_obtained_beforedetection];
pfacolbeforedetection3=[pfacolbeforedetection3 pfa_obtained_beforedetection];
pdcolafterdetection3=[pdcolafterdetection3 pd_obtained_afterdetection];
pfacolafterdetection3=[pfacolafterdetection3 pfa_obtained_afterdetection];
end
figure
plot(pfacolbeforedetection3,pdcolbeforedetection3,’*’)
figure
plot(pfacolafterdetection3,pdcolafterdetection3,’*’)
pdcolbeforedetection4=[];
pfacolbeforedetection4=[];
pdcolafterdetection4=[];
pfacolafterdetection4=[];
x=[0.4 0.2];
for p0=0.01:1/1000:0.99
[datatx,datarx,datadetected,bayesruleselected,pl,ph,val_pl,val_ph,...
pd,gamma,pd_obtained_beforedetection,pfa_obtained_beforedetection,...
pd_obtained_afterdetection,pfa_obtained_afterdetection]...
=binarychanneldetection(x,y,z,len,p0,pfa);
pdcolbeforedetection4=[pdcolbeforedetection4 pd_obtained_beforedetection];
pfacolbeforedetection4=[pfacolbeforedetection4 pfa_obtained_beforedetection];
pdcolafterdetection4=[pdcolafterdetection4 pd_obtained_afterdetection];
pfacolafterdetection4=[pfacolafterdetection4 pfa_obtained_afterdetection];
end
figure
plot(pfacolbeforedetection4,pdcolbeforedetection4,’*’)
figure
plot(pfacolafterdetection4,pdcolafterdetection4,'*')