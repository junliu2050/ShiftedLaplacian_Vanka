%Solve 2D Shifted Poisson equation with Geometric Multgrid method
clear
set(0, 'defaultaxesfontsize',20,'defaultaxeslinewidth',1,...
    'defaultlinelinewidth',2,'defaultpatchlinewidth',1.5,...
    'defaulttextfontsize',20,'defaulttextInterpreter','latex');
global LEVEL0;
fprintf('N\t Rel. Res.\t  Error \t iter \t Rate\n');
T=1; maxit=100;tol=1e-10; mu1=1;mu2=0;  LEVEL0=3; %3
smlist={'Jacobi','Jacobi-C','Vanka','Mass','Mass-shift','Vanka-e','Vanka-e2','GS','RB-GS'};
lsty ={'x-','d-','^-','v-','>-','<-','p-','h-','.-'};
for kk=[1 3] %:length(smlist)
    smtype=smlist{kk}; %select compared smoothers
    fprintf('--------------Smoother type: [%s]-------------------\n',smtype)
    for n=6:6
        nt=2^(n-1);
        %uncomment to select different examples
        %B from Direct PinT solver (heat equation case)
        exname='Heat';  dt=T/nt; eigB=eigBNewton(nt)/dt; %Example 1
        %exname='Wave'; dt=T/nt; eigB=(eigBNewton(nt)/dt).^2; %not working
        
        %B from Backward Heat Conduction Problem paper (PinT-QBVM) 
        %     delta=1e-2;alpha=delta; exname='BackwardHeat'; %Example 2
        %     e=ones(nt+1,1);
        %     B=spdiags([-e e]/dt,[-1 0],nt+1,nt+1);
        %     B(1,end)=1/alpha; eigB=eig(full(B));
        %exname='Helmholtz'; eigB=-(0:nt).^2*(1-0.5i); %Example 3
        tic
        iterlist=[]; ratelist=[];
        for jj=1:nt
            lamj=eigB(jj);
            fineL=n;
            xmin=0;xmax=1;
            ymin=0;ymax=1;
            %Setup operators
            mg=[];
            %define the matrix at each level
            for Level=fineL:-1:LEVEL0
                nn=2^Level; h=(xmax-xmin)/nn;m=nn-1;
                [mg(Level).A]=(1/h^2)*gallery('poisson',m)+(lamj)*speye(m^2); %operator in all levels
                e=ones(m,1);
                Me=h*spdiags([e 4*e e]/6,-1:1,m,m); % 1D mass matrix
                %interpolation operator in matrix form
                Pn=(1/2)*spdiags([e 2*e e],-2:0,m,m);
                mg(Level).P=kron(Pn(:,1:2:end-2),Pn(:,1:2:end-2));

                mg(Level).M=[]; E=speye(m^2); mg(Level).shift=0;
                switch(smtype)
                    case 'Jacobi-C'
                        %Thm 4.2 in https://www.cs.ubc.ca/~greif/Publications/hg2021.pdf
                        c=2/h^2;cm=1/h^2;Lam=4/h^2+lamj;
                        b1=1-2*(c-cm)/Lam;b2=1+2*c/Lam;
                        mg(Level).shift=(abs(b1)/b1+abs(b2)/b2)/(abs(b1)+abs(b2));
                        mg(Level).factor=abs(b1-b2)/(abs(b1)+abs(b2));
                    case 'Vanka'
                        %mg(Level).M=@(r) Vanka2D_e(r,h,lamj); % matrix free M*r
                        %matrix version for faster CPU times
                        z=h^2*lamj;%complex number
                        aa=(z^2 + 8*z + 14)/(z^3 + 12*z^2 + 44*z + 48);
                        bb=1/(z^2 + 8*z + 12);cc=2/(z^3 + 12*z^2 + 44*z + 48);
                        Me=spdiags(h*[sqrt(cc)*e (2*bb/sqrt(cc))*e sqrt(cc)*e]/2,-1:1,m,m);
                        %mg(Level).M=kron(Me,Me)+h^2*(aa-bb^2/cc)*E;%large matrix
                        mg(Level).M=Me; mg(Level).shift=h^2*(aa-bb^2/cc);%small matrix
                    case 'Mass-shift'
                        z=h^2*lamj;
                        aa=(z + 2)/(z^2 + 4*z + 3);bb=1/(z^2 + 4*z + 3);
                        Me=h*spdiags([bb*e 2*aa*e bb*e]/2,-1:1,m,m); %shifted 1D mass matrix
                        mg(Level).M=kron(Me,Me);
                    case 'Vanka-e'
                        mg(Level).M=(3*kron(Me,Me)+h^2*E)/8; %Me in (23)
                    case 'Vanka-e2'
                        % Me2=h*spdiags([e 3*e e],-1:1,m,m); %
                        %mg(Level).M=(kron(Me2,Me2)+h^2*5*E)/28;
                        Me3=h*spdiags([e (10/3)*e e],-1:1,m,m); %optimal
                        mg(Level).M=(3*kron(Me3,Me3)+h^2*(32/3)*E)/80;
                    case 'Mass'
                        mg(Level).M=kron(Me,Me); %Mfe in (22) s
                end
            end

            A=mg(fineL).A;
            %Setup right-hand-side
            N=2^fineL;
            hx=(xmax-xmin)/N;
            hy=(ymax-ymin)/N;
            xgrid=xmin+(1:N-1)'*hx;%interior nodes
            ygrid=ymin+(1:N-1)'*hy;
            [X,Y] = meshgrid(xgrid,ygrid);
            S0=U0(X,Y);%true solution exact
            %For original equaion
            F=fF(X,Y);
            %F(1,:)=0;F(end,:)=0;F(:,1)=0;F(:,end)=0;
            b=F(:);
            %Begin Multigrid Algorithm
            x=rand(size(b));
            h=hx;
            r0=h*norm(b-A*x);rk=r0;
            err0=[]; it=1; err0(it)=rk;
            %x=full_mg_2d(mg,x,b,fineL,1,1,smtype);%Full Multigrid, faster
            while((it<=maxit)&&(rk>=r0*tol))
                it=it+1;
                x=mg_iter_2d(mg,x,b,fineL,mu1,mu2,smtype);                
                rk=h*norm(b-A*x);
                err0(it)=rk;
            end
            iter=length(err0)-1;%total V-cycle iterations
            %use AA-MG
            %[x,iter,err0]=AndAcc(@(x)mg_iter_2d(mg,x,b,fineL,mu1,mu2,smtype),x,10,maxit,tol,tol);
            
            res=err0(end)/r0;
            rate=res^(1/iter);
            %semilogy(0:length(err0)-1,err0/r0,lsty{kk},'DisplayName',sprintf('%s(CPU=%1.2f s, iter=%d, rate=%1.3f)',smtype,cpu,iter,rate));axis tight; hold on
            %plot(jj,iter,lsty{kk}); hold on
            iterlist=[iterlist;iter]; ratelist=[ratelist;rate];
            err=norm(S0(:)-x,inf);%solution max error, decrease by 4 times

            fprintf('[%d/%d]&\t%d&\t %1.1e & \t %1.1e &\t %d& \t%1.3f \\\\ \n',jj,nt,N,res,err,iter,rate);
        end
        cpu=toc;%CPU time should grows linearly, 4 times
        subplot(1,2,1);axis tight
        plot((1:nt),iterlist,lsty{kk},'DisplayName',sprintf('%s(CPU=%1.2f s)',smtype,cpu)); hold on
        xlabel(sprintf('system index ($1\\le j\\le %d$)',nt));title('Number of Iterations');
        subplot(1,2,2);axis tight
        plot((1:nt),ratelist,lsty{kk},'DisplayName',sprintf('%s(max rate=%1.3f)',smtype,max(ratelist))); hold on
        xlabel(sprintf('system index ($1\\le j\\le %d$)',nt));title('Convergence rate');
    end
end
set(gcf, 'Position', get(0, 'Screensize').*[1 1 2/3 1/2]);
subplot(1,2,1);legend('-DynamicLegend','Location', 'Best');
subplot(1,2,2);legend('-DynamicLegend','Location', 'Best');
%sgtitle(sprintf('%d-level GMG with W(%d+%d) cycles $(h=1/%d,\\tau=1/%d)$',(fineL-LEVEL0)+1,mu1,mu2,N,nt));
%xlabel('V-cycle Iterations'); ylabel('Relative residual norms')
figname1=sprintf('../Manuscript/%s_mesh_%d_T%d',exname,nt,100*T);
%print(figname1,'-depsc')
%  figure(2)
% plot(eigB,'+');
% set(gcf, 'Position', get(0, 'Screensize').*[1 1 2/3 1/2]);
% %title(sprintf('All [%d] eigenvalues of $B$',length(eigB)));
% title('The distribution of all complex shifts'); xlabel('Real parts');ylabel('Imaginary parts');
% figname2=sprintf('../Manuscript/%s_eigB_%d_T%d',exname,nt,100*T);
% print(figname2,'-depsc')


function [x]=mg_iter_2d(mg,x0,b,level,pre,post,smoother)
%fprintf('+++++++++++++++++++Level[%d] Begin+++++++++++++++++++++\n',level);
A=mg(level).A;
global LEVEL0;
if(level==LEVEL0) %Coarest level
    x=A\b;
else
    % presmooth 
    x = mg_smooth(A,x0,b,pre,smoother,mg(level).M,mg(level).shift); 
    % Restrict residual
    r = b - A*x;
    %rc = mg_restricth(r,h);
    %rc = rest(r);
    rc=(mg(level).P)'*r/4;
    % coarse grid correction
    cc = mg_iter_2d(mg,zeros(size(rc)),rc,level-1,pre,post,smoother);
    cc = mg_iter_2d(mg,cc,rc,level-1,pre,post,smoother); %add this for W cycle
    %x = x + intp2(cc);
    %x=x + mg_interph(cc,h/2);
    x = x + mg(level).P*cc; %
    % postsmooth
    %fprintf('Level[%d],Before Postsmooth res=%1.2e\n',level,norm(b-A*x));
    x = mg_smooth(A,x,b,post,smoother,mg(level).M,mg(level).shift);
    %fprintf('Level[%d],After Postsmooth res=%1.2e\n',level,norm(b-A*x));
end
%fprintf('+++++++++++++++++++Level[%d] End+++++++++++++++++++++\n',level);
end
function [x]=mg_smooth(A,x0,b,nv,smoother,M,shift)
x=x0; 
switch smoother
    case {'Mass','Mass-shift'} %mass smoother
        w=3/4;
        for k = 1:nv
            x=x + w*M*(b - A*x);
        end
    case 'Vanka'
        %w=24/25;
        r=shift; w=(20*r^2+20*r+48)/(2*r^2+21*r+50);
        for k = 1:nv
            res=(b - A*x);
            %x=x + w*M*(b - A*x); %avoid construct large M matrix
            Me_res=M*reshape(res,sqrt(length(b)),[])*M;
            x=x+w*(Me_res(:)+shift*res);
        end
    case 'Vanka-e'
        %w=24/25;
        r=shift; w=(20*r^2+20*r+48)/(2*r^2+21*r+50);
        for k = 1:nv
            x=x + w*M*(b - A*x);
        end
    case 'Vanka-e2'
        w=0.5253; %w=0.5801;
        for k = 1:nv
            x=x + w*M*(b - A*x);
        end
    case 'Jacobi'
        w=4/5;
        dinv = 1./diag(A);
        for k = 1:nv
            x=x + w*dinv.*(b - A*x);
        end
    case 'Jacobi-C'
        w=shift;
        dinv = 1./diag(A);
        for k = 1:nv
            x=x + w*dinv.*(b - A*x);
        end
    case 'GS' %G-S Smoother
        L = tril(A);
        for k = 1:nv
            x=x + (L\(b - A*x));
        end
    case 'RB-GS'
        dinv = 1./diag(A);
        for k = 1:nv
            x(1:2:end)=x(1:2:end) + dinv(1:2:end).*(b(1:2:end) - A(1:2:end,:)*x);
            x(2:2:end)=x(2:2:end) + dinv(2:2:end).*(b(2:2:end) - A(2:2:end,:)*x);
        end
end
end 

function z=fF(x,y)
%z=2*pi^2*sin(pi*x).*sin(pi*y);
z=2*(1-6*x.^2).*(y.^2-y.^4)+2*(1-6*y.^2).*(x.^2-x.^4);
end
function z=U0(x,y)
%z=sin(pi*x).*sin(pi*y);
z=(x.^2-x.^4).*(y.^4-y.^2);
end