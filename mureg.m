
function [coe, res]=mureg(xx, yy)

coe=(xx'*xx)^-1*xx'*yy;

zz=xx*coe-yy;

res=sqrt(mean(zz.^2))/sqrt(mean(yy.^2));

figure;
plot(yy,'r','LineWidth', 2);
hold on
plot(xx*coe, 'b', 'LineWidth',1);
legend('Red-Real Data ', 'Blue-Regession Data');
px=double(int8(length(yy)*2/3));
py=double(abs(int8((max(yy)-min(yy))/2)));
text(px,py,['Residual: ' num2str(res)],'HorizontalAlignment','left')

display(['Residual: ', num2str(res)]);

