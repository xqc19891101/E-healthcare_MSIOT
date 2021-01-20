function [number,path] = findPath(mark,start,stop)
if length(mark{stop,1}) == 1
    number = 1;
    if mark{stop,1} ~= start
         T=mark{stop,1};
         [number,path1] = findPath(mark,start,T);
          for n = 1:number
              path{n} =[ path1{n},stop];
          end
    else
        path = {[start,stop]};
    end
else 
    sNo = length(mark{stop,1});
    number = sNo;
    path = cell(sNo,1);
    for n = 1:sNo
         T = mark{stop,1}(n);
         [number,path1] = findPath(mark,start,T);
          for k = 1:number
              path{n} =[ path1{k},stop];
          end
    end    
end

end