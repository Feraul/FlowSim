function [vec_gravelem,vec_gravface,vec_gravpoint,gravelem,gravpoint,gravface]=PLUG_Gfunction
%function [grav]=PLUG_Gfunction
global  numcase g elem centelem centface bedge inedge coord u0

vec_gravelem=zeros(size(elem,1),3);
vec_gravface=zeros(size(bedge,1)+size(inedge,1),3);
vec_gravpoint=zeros(size(coord,1),1);
gravelem=zeros(size(elem,1),1);
gravpoint=zeros(size(coord,1),1);
gravface=zeros(size(bedge,1)+size(inedge,1),1);
switch numcase
    
    case 18
        % calculo da gravidade no baricentro
        for i=1:size(elem,1)
            y = centelem(i,2);
            gravelem(i,:)=g(2,1)*y;
            
        end
        % calculo da gravidade nas interfaces
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
                
                a=0.5*(coord(v1,:)+coord(v2,:));
                y2=a(1,2);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
                
                a=0.5*(coord(v1,:)+coord(v2,:));
                y2=a(1,2);
            end
            gravface(j)= g(2,1)*y2 ;
        end
        
        % 4.2 Numerical Examples: incompressible flow in a unit square domain.
        % The domain has a discontinuity line and the gravitational forces are
        % given as a linear combination of two contributions
        % STARNONI et al. Consistent MPFA Discretization for Flow in the
        % Presence of Gravity. 2019.
        
        % 4.2 Test 1 where a1=1 e a2=0
    case 421
        
        % Para r=0;
        delta=0.5;
        
        h1=10;
        h2=1;
        
        for i=1: size(elem,1)
            
            y = centelem(i,2);
            
            
            if(y>delta)
                % gravity vector
                grav(i,1)=0;
                grav(i,2)=h1;
                
                % gravity through the edge element
                gravelem(i)= h1*y - cte;
                
            else
                % gravity vector
                grav(i,1)=0;
                grav(i,2)=h2;
                
                % gravity through the edge element
                gravelem(i)= h2*y - cte;
            end
            
        end
        
        for j=1: size(coord, 1)
            
            %Define "x" and "y"
            y1 = coord(j,2);
            
            if(y1>delta)
                gravpoint(j,1)= h1*y1 - cte;
            else
                gravpoint(j,1)= h2*y1 - cte;
            end
            
        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
                
                a=0.5*(coord(v1,:)+coord(v2,:));
                y2=a(1,2);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
                
                a=0.5*(coord(v1,:)+coord(v2,:));
                y2=a(1,2);
            end
            
            if(y2>delta)
                gravface(j)= h1*y2 - cte;
            else
                gravface(j)= h2*y2 - cte;
            end
            
        end
        
        
        
        
        % 4.2 Test 2 where a1=0 e a2=1
    case 422
        
        for i=1: size(elem,1)
            
            x = centelem(i,1);
            y = centelem(i,2);
            
            % gravity vector
            vec_gravelem(i,1)=cos(x)*cos(y);
            vec_gravelem(i,2)=-sin(x)*sin(y);
            
            
            % gravity through the edge element
            
            gravelem(i,1)=-sin(x)*cos(y)-u0;
            
            
        end
        
        %gravpoint
        for j=1:size(coord,1)
            %Define "x" and "y"
            x1 = coord(j,1);
            y1 = coord(j,2);
            % parametro segundo  Starnoni
            
            % gravity vector
            vec_gravpoint(j,1)=cos(x1)*cos(y1);
            vec_gravpoint(j,2)=-sin(x1)*sin(y1);
            
            % solucao analitica
            gravpoint(j,1)= -sin(x1)*cos(y1)-u0;
        end
        
        %gravface
        for j=1:size(bedge,1)+size(inedge,1)
            
            x2=centface(j,1);
            y2=centface(j,2);
            %Define "x" and "y"
            %             if j<=size(bedge,1)
            %                 v1=bedge(j,1);
            %                 v2=bedge(j,2);
            %
            %                 a=0.5*(coord(v1,:)+coord(v2,:));
            %                 x2=a(1,1);
            %                 y2=a(1,2);
            %             else
            %                 v1=inedge(j-size(bedge,1),1);
            %                 v2=inedge(j-size(bedge,1),2);
            %
            %                 a=0.5*(coord(v1,:)+coord(v2,:));
            %                 x2=a(1,1);
            %                 y2=a(1,2);
            %             end
            vec_gravface(j,1)=cos(x2)*cos(y2);
            vec_gravface(j,2)=-sin(x2)*sin(y2);
            
            gravface(j,1)= -sin(x2)*cos(y2)-u0;
            
        end
        
   
        % 4.2 Test 3 where a1=1 e a2=1
    case 423
        
        % Para r=0;
        delta=0.5;
        
        h1=10;
        h2=1;
        
        for i=1: size(elem,1)
            
            x = centelem(i,1);
            y = centelem(i,2);
            
            if(y>delta)
                % gravity vector
                grav(i,1)=-cos(x)*cos(y);
                grav(i,2)=h1+sin(x)*sin(y);
                
                % gravity through the edge element
                gravelem(i)= - h1*y + sin(x)*cos(y) - cte;
                
            else
                % gravity vector
                grav(i,1)=-cos(x)*cos(y);
                grav(i,2)=h2+sin(x)*sin(y);
                
                % gravity through the edge element
                gravelem(i)= - h2*y + sin(x)*cos(y) - cte;
                
            end
            
        end
        
        for j=1: size(coord, 1)
            
            %Define "x" and "y"
            x1 = coord(j,1);
            y1 = coord(j,2);
            % parametro segundo  Starnoni
            
            % solucao analitica
            if(y1>delta)
                gravpoint(j)= - h1*y1 + sin(x1)*cos(y1) - cte;
            else
                gravpoint(j)= - h2*y1 + sin(x1)*cos(y1) - cte;
            end
            
        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
                
                a=0.5*(coord(v1,:)+coord(v2,:));
                x2=a(1,1);
                y2=a(1,2);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
                
                a=0.5*(coord(v1,:)+coord(v2,:));
                x2=a(1,1);
                y2=a(1,2);
            end
            
            if(y2>delta)
                gravface(j)= -h1*y2 + sin(x2)*cos(y2) - cte;
            else
                gravface(j)= -h2*y2 + sin(x2)*cos(y2) - cte;
            end
        end
        
        
        % 4.2 Test 4 where a1=1 e a2=100
    case 424
        
        % Para r=0;
        delta=0.5;
        
        h1=10;
        h2=1;
        
        for i=1: size(elem,1)
            
            x = centelem(i,1);
            y = centelem(i,2);
            
            if(y>delta)
                grav(i,1)= -100*cos(x)*cos(y);
                grav(i,2)=h1+100*sin(x)*sin(y);
                
                gravelem(i,1)=100*sin(x)*cos(y)-h1*y;
            else
                grav(i,1)=-100*cos(x)*cos(y);
                grav(i,2)=h2+100*sin(x)*sin(y);
                
                gravelem(i,1)=100*sin(x)*cos(x)-h2*y;
            end
            
        end
        
        for j=1: size(coord, 1)
            
            %Define "x" and "y"
            x1 = coord(j,1);
            y1 = coord(j,2);
            % parametro segundo  Starnoni
            
            % solucao analitica
            if(y1>delta)
                gravpoint(j,1)= -h1*y1 + 100*sin(x1)*cos(y1) - cte;
            else
                gravpoint(j,1)= -h2*y1 + 100*sin(x1)*cos(y1) - cte;
            end
            
        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
                
                a=0.5*(coord(v1,:)+coord(v2,:));
                x2=a(1,1);
                y2=a(1,2);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
                
                a=0.5*(coord(v1,:)+coord(v2,:));
                x2=a(1,1);
                y2=a(1,2);
            end
            
            if(y2>delta)
                gravface(j)= -h1*y2 + 100*sin(x2)*cos(y2) - cte;
            else
                gravface(j)= -h2*y2 + 100*sin(x2)*cos(y2) - cte;
            end
        end
    case 425
        for i=1: size(elem,1)
            grav(i,1)= g(1,1);
            grav(i,2)= g(2,1);
            grav(i,3)= g(3,1);
        end
        
    otherwise
        %==================================================================
        % quando a gravidade constante
        for i = 1:size(centelem,1)
            % calculo da gravidade no elemento
            vec_gravelem(i,:)=g';
        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            % calculo da gravidade na interface
            vec_gravface(j,:)= g';
        end
        
end


