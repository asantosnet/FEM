close all

% We will calculate using the three methods and plot the Convergence Curve

TOL = 10e-6;

[ test1xlastBFSG,test1ylastBFSG,test1xBFSG,test1yBFSG,test1xelasticBFSG,test1yelasticBFSG ] = HP_BFGS( );

[ test1xlastN,test1ylastN,test1xN,test1yN,test1xelasticN,test1yelasticN ] = HP_Newton( );

[ test1xlastMN,test1ylastMN,test1xMN,test1yMN,test1xelasticMN,test1yelasticMN ] = HP_Modified_Newton(  );


figure (11)


% We will plot using log scale

loglog(test1xBFSG,test1yBFSG,'-b');


% We will force both axis to be proportional, which is not achieveable by
% using equal axis when ploting with loglog

xLimits = [TOL*10e-3 10e5];
yLimits = [TOL*10e-3 10e5];
axes_equal_loglog (gca,xLimits,yLimits)

xlabel('|g(k)|');
ylabel('|g(k+1)|');

hold on

loglog(test1xlastBFSG,test1ylastBFSG,'-r');
loglog(test1xelasticBFSG,test1yelasticBFSG,'-g');

loglog(test1xN,test1yN,'--b');
loglog(test1xelasticN,test1yelasticN,'--g');
loglog(test1xlastN,test1ylastN,'--r');

loglog(test1xMN,test1yMN,'-.b');
loglog(test1xlastMN,test1ylastMN,'-.r');
loglog(test1xelasticMN,test1yelasticMN,'-.g');

grid on

legend(' BFSG  ','BFSG For last step  ','BFSG Elastic domain',...
    'Newton ','Newton  For last step ','Newton  Elastic domain',...
    ' MNewton ','MNewton  For last step ','MNewton  Elastic domain');

hold off