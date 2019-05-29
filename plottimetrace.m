load stab86166.dat
subplot(2,1,1);                         
plot(stab86166(:,1),stab86166(:,3),'b');    
hold on                                 
plot(stab86166(:,1),stab86166(:,4),'r');
plot(stab86166(:,1),stab86166(:,5),'g');
plot(stab86166(:,1),stab86166(:,6),'k');
text(2900,-0.28,'D_I','Color','b');
text(2900,0.2,'D_R','Color','r');
text(2900,0.8,'D_{nc}','Color','g');
text(2900,0.0,'D_{tot}','Color','k');
hold off
subplot(2,1,2);
plot(stab86166(:,1),stab86166(:,7),'k');
hold on                                 
plot(stab86166(:,1),stab86166(:,8),'r');
xlabel('Time (ms)');
text(2900,2.8,'\beta_N','Color','k');
text(2900,3.7,'4*li','Color','r');   
hold off

