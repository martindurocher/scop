function fit = estimate(obj,method,data)

%default parameters
if(nargin < 2); method = 'pairwise'; end
       
switch method
    case 'ind'
        fit = obj;
    case 'pairwise'
        if(nargin < 3); 
            data = obj.pairs; 
        end
                
        pid = obj.pdistance < obj.HTOL;
        u0 = data(pid,:);
        h0 = obj.pdistance(pid);    
        fit = obj.sp.fit('pairwise',h0,u0);
            
    case 'full'
        if(nargin < 3); 
            data = obj.ranks; 
        end
            
        fit = obj.sp.fit('full',obj.distance,data);
           
end
