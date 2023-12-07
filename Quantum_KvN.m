clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%  Parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%

nx=2^5 ; %nuber of grid points
init_d=3; %Diplacement of the original state

%grid parameters
Lx=(nx-1)*(2*pi/nx)^.5; Lx0=-Lx/2; dx=Lx/(nx-1);x=Lx0+dx*(0:nx-1);
dpx=(2*pi)/(nx*dx); px=dpx*(0:nx/2); px=[px(1:end-1),-px(end:-1:2)];

%freqs and couplings
omega=10; omega_a=4.403989664625457e-05;omega_b=48*omega_a;
ga=omega; gb=omega*1.4;


%SIMULATION TIME PARAMETERS 
number_of_oscillations=1.5; trotter_steps=300;
tf=number_of_oscillations*2*pi/omega
dt=tf/trotter_steps; 
n=round(tf/dt);  
tt=linspace(0,tf,n);
dt=tt(2)-tt(1); 


%Boolean variables
plot_dyamics=0; %Plot the dynamics
track_quantities=1; %Track Purities and LN (is quite slow for 7 qubits, 3 days)

P_spin_a=[]; P_spin_b=[]; P_spin_ho=[]; P_spins=[]; %PURITIES
LNs=[]; %LOGARITHMIC NEGATIVITY 

 

%%%%%%%%%%%%%%%%%%%%%% INITIAL STATE %%%%%%%%%%%%%%%%%%%%%%

%Note: The encoding is a vector of matrices to identify  the subspaces better

%x and p  
psi0_q=exp(-(x-Lx/init_d).^2/2)/norm(exp(-(x-Lx/init_d).^2/2));
psi0_p=fft(exp(-(x-Lx/init_d).^2/2))/norm(fft(exp(-(x-Lx/init_d).^2/2)));


% figure
% plot(px,abs(psi0_p))
psi0_p=[psi0_p(nx/2+1:end), psi0_p(1:nx/2)];

psi0_qp=psi0_q'*psi0_p; %(matrix encoding)  

%spins
psi0_spin_a=1/2^.5*[1,1i];
psi0_spin_b=1/2^.5*[1,1i];

%total
psi0=kron(psi0_spin_a,kron(psi0_spin_b,psi0_qp));

psi=psi0;



%%%%%%%%%%%%%%%%%%  OPERATORS DEFINITION %%%%%%%%%%%%%%%%% 

for k=1:nx
    x_p(k,:)=[px(nx/2+1:end), px(1:nx/2)]; 
    l_p(k,:)=[x(nx/2+1:end), x(1:nx/2)];
    l_q(:,k)=px';
    x_q(:,k)=x';
end


lq_p=l_q.*x_p;
q_lp=-l_p.*x_q;

%harmonic oscillator
kin=omega*kron([1,1,1,1], lq_p);
expkin=exp(-1i*dt*kin);
V=omega*kron([1,1,1,1], q_lp); 
expV=exp(-1i*V*dt);

%spins
spin_a=omega_a*[ones(nx,2*nx),-ones(nx,2*nx)];
spin_b=omega_b*[ones(nx,nx), -ones(nx,nx), ones(nx,nx), -ones(nx,nx)];

%coupling
spin_a_lp=(2)^.5*ga*kron([1,-1],kron([1,1],l_p));
spin_b_lp=(2)^.5*gb*kron([1,1],kron([1,-1],l_p));


  
%%%%%%%%%%%%%% Initial track of the quantities  %%%%%%%%%%%%%% 

if track_quantities==1

     psi_vect=[]; %we must reshape psi into a vector

     for  k=1:4
        psi_aux=psi(:,(k-1)*nx+1:k*nx);

        for f=1:nx
             psi_vect=[psi_vect, psi_aux(f,:)];
        end 
    end 

    psi_vect=psi_vect.';
    rho= psi_vect*psi_vect';

    %Purity
    XPT_spin_a = PartialTrace(rho,[2,3],[2,2,nx*nx],-1);
    XPT_spin_b = PartialTrace(rho,[1,3],[2,2,nx*nx],-1);
    XPT_spins = PartialTrace(rho,[2],[4,nx*nx],-1);
    XPT_ho = PartialTrace(rho,[1,2],[2,2,nx*nx],-1);
    P_spin_a(1)=trace(XPT_spin_a* XPT_spin_a);
    P_spin_b(1)=trace(XPT_spin_b* XPT_spin_b);
    P_spin_ho(1)=trace(XPT_ho* XPT_ho);
    P_spins(1)=trace(XPT_spins* XPT_spins);

    %Logarithmic Negativity between spins
    rho_12 = PartialTrace(rho,2,[4,nx*nx],-1);
    rho_T2 = PartialTranspose(rho_12,2,[2,2]);
    LN=log2(sum(eig(rho_T2*rho_T2').^.5));
    LN=max([real(LN),0]);
    LNs(1)=LN;
    
%         %Logarithmic Negativity between spin and qho
%         rho_reduced = PartialTrace(rho,[2],[2,2,nx*nx],-1);
%         rho_T2 = PartialTranspose(rho_reduced,2,[2,nx*nx]);
%         LN=log2(sum(eig(rho_T2*rho_T2').^.5));
%         LN=max([real(LN),0]);
%         LNs_qs(1)=LN;

end
    



%%%%%%%%%%%%%%%% Time evolution %%%%%%%%%%%%%%%%%
    
    for i=1:n-1

        disp(i)

        psi=psi.*exp(-1i*spin_a*dt).*exp(-1i*spin_b*dt); 
        psi_lp=[];
        
        for f=1:4
             psi_lp=[psi_lp,fft(psi(:,(f-1)*nx+1:f*nx), nx, 2)];
        end 
       
    
        psi_lp=psi_lp.*expV.*exp(-1i*(spin_a_lp+spin_b_lp)*dt);
    
        
        psi=[];
    
        for f=1:4
             psi=[psi,ifft(psi_lp(:,(f-1)*nx+1:f*nx), nx, 2)];
        end 
    
        psi_lq=[];
    

        for f=1:4
             psi_lq=[psi_lq,fft(psi(:,(f-1)*nx+1:f*nx), nx, 1)];
        end 
    
    
        psi_lq=psi_lq.*expkin;
        
        
        psi=[];
    
        for f=1:4
             psi=[psi,ifft(psi_lq(:,(f-1)*nx+1:f*nx), nx, 1)];
        end 
    

        psi_plot=zeros(nx);
        
        
        if plot_dyamics==1

            for f=1:4
               psi_plot=psi_plot+abs(psi(:,(f-1)*nx+1:f*nx)).^2;
            end 
            imagesc(abs(psi_plot))
            set(gca,'YDir','normal')
            drawnow

        end
    
    
        
        if track_quantities==1

            psi_vect=[];
        
            for  k=1:4
        
              psi_aux=psi(:,(k-1)*nx+1:k*nx);
        
              for f=1:nx
        
                 psi_vect=[psi_vect, psi_aux(f,:)];
        
              end 
        
            end 
            
            psi_vect=psi_vect.';
            rho= psi_vect*psi_vect';

            %Spins purity
            XPT_spin_a = PartialTrace(rho,[2,3],[2,2,nx*nx],-1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
            XPT_spin_b = PartialTrace(rho,[1,3],[2,2,nx*nx],-1);
            XPT_ho = PartialTrace(rho,[1,2],[2,2,nx*nx],-1);
            XPT_spins = PartialTrace(rho,[2],[4,nx*nx],-1);
            P_spin_a(i+1)=trace(XPT_spin_a* XPT_spin_a);
            P_spin_b(i+1)=trace(XPT_spin_b* XPT_spin_b);
            P_spin_ho(i+1)=trace(XPT_ho* XPT_ho);  
            P_spins(i+1)=trace(XPT_spins* XPT_spins);   
             
            %Ln betwen spins
            rho_12 = PartialTrace(rho,2,[4,nx*nx],-1);
            rho_T2 = PartialTranspose(rho_12,2,[2,2]);
            LN=log2(sum(eig(rho_T2*rho_T2').^.5));
            LN=max([real(LN),0]);
            LNs(i+1)=LN;
        
    
        %Ln between spin a and ho
%         rho_reduced = PartialTrace(rho,[2],[2,2,nx*nx],-1);
%         rho_T2 = PartialTranspose(rho_reduced,2,[2,nx*nx]);
%         LN=log2(sum(eig(rho_T2*rho_T2').^.5));
%         LN=max([real(LN),0]);
%         LNs_qs(1+i)=LN;
       end 
    
    end

    if track_quantities==1

        colors=[[0.97,0.91,0.99];[0.85,0.76,0.89];[0.67,0.53,0.73];[0.51,0.33,0.59];[0.6836,    0.6641,    0.6914]; [0.44,0.23,0.53];[0.30,0.09,0.39]];
        figure(2)
        hold on
        
        plot(linspace(0,1.5,trotter_steps),P_spin_a, '--','color',colors(7,:),'LineWidth',1.5)
        plot(linspace(0,1.5,trotter_steps),P_spin_b, '-','color',colors(5,:),'LineWidth',1.5)
        plot(linspace(0,1.5,trotter_steps),P_spin_ho, '-.','color',colors(3,:),'LineWidth',1.5)
        legend({'Purity_spin_a','Purity_spin_b','Purity_ho'},'Interpreter','latex','FontSize',10)
        xlabel('$t/\omega$', 'FontSize',20,'Interpreter','latex')
        ylabel('Purity', 'FontSize',20,'Interpreter','latex','Color','k')
        legend boxoff
        
        figure(3)
        plot(linspace(0,1.5,trotter_steps),LNs,'-','color',colors(3,:),'LineWidth',1.5)
        xlabel('$t/\omega$', 'FontSize',20,'Interpreter','latex')
        ylabel('Logarithmic Negativity', 'FontSize',20,'Interpreter','latex','Color','k')
        legend({'Logarithmic negativity'},'Interpreter','latex','FontSize',20)
        hold on
        ax = gca;
        ylim([-0.1,0.7])
        legend boxoff
        
    end
   