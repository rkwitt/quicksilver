M = 8;
K = 50;
h = figure;
p1 = plot(1:1:K*M,avar(1:K*M)/avar(1),'bo-','MarkerFaceColor','b');
hold on;
p2 = plot(1:M:K*M,avar(1:M:K*M)/avar(1),'r*:', 'MarkerSize',10);
title({'LDMM3D Convergence';sprintf('%d Spheres; %d Iterations',M,K)});
xlabel('Sub-Iteration');
ylabel('Integrated Pointwise Variance (Normalized)');
legend('Sub-Iteration','Iteration');

boldifyPlot(h,p1);
boldifyPlot(h,p2);