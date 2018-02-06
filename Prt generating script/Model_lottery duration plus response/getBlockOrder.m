function block_order = getBlockOrder(is_mon, Datamon, Datamed, tpb)
%BLOCK_ORDER Determine Condition Order
% Which blocks contain the required domain, across the entirety of the
% experiment? This is required because the imaging file isn't segregated
% between mon and med. This function extracts them from collected data.
%
%OUTPUT 
% block_order - 1x4 cell that contains the IDs of blocks in a given domain.

% Gets 1x3 matrix of hour, minute, second, of the first trial of each
% domain

%% for design of alternating between mon and med in two blocks
monstart = Datamon.trialTime(1).trialStartTime(4:6);
medstart = Datamed.trialTime(1).trialStartTime(4:6);

if is_mon % domain is monetary
    if monstart(1)*60 + monstart(2) > medstart(1)*60 + medstart(2) %time in minutes. if medical is first
        block_order{1} = '3'; % medical first
        block_order{2} = '4';
        block_order{3} = '7';
        block_order{4} = '8';
    else
        block_order{1} = '1'; % monetary first
        block_order{2} = '2';
        block_order{3} = '5';
        block_order{4} = '6';       
    end
else % domain is medical
    if monstart(1)*60 + monstart(2) > medstart(1)*60 + medstart(2) % if medical is first
        block_order{1} = '1'; % medical first
        block_order{2} = '2';
        block_order{3} = '5';
        block_order{4} = '6';        
    else
        block_order {1}= '3'; % monetary first
        block_order{2} = '4';
        block_order{3} = '7';
        block_order{4} = '8';

    end
end


%% below are the scripts used for random block assigment.
% if is_mon % domain is monetary
%     if monstart(1)*60 + monstart(2) > medstart(1)*60 + medstart(2) %time in minutes. if medical is first
%         block_order{1} = '2'; % medical first
%     else
%         block_order{1} = '1'; % monetary first
%     end
% else % domain is medical
%     if monstart(1)*60 + monstart(2) > medstart(1)*60 + medstart(2)
%         block_order{1} = '1'; % medical first
%     else
%         block_order {1}= '2'; % monetary first
%     end
% end
% 
% % NOTE on index: tpb+1= the first trial of the second block.
% monstart = Datamon.trialTime(tpb+1).trialStartTime(4:6);
% medstart = Datamed.trialTime(tpb+1).trialStartTime(4:6);
% 
% if is_mon % domain is monetary
%     if monstart(1)*60 + monstart(2) > medstart(1)*60 + medstart(2) %time in minutes
%         block_order{2} = '4'; % medical first
%     else
%         block_order{2} = '3'; % monetary first
%     end
% else % domain is medical
%     if monstart(1)*60 + monstart(2) > medstart(1)*60 + medstart(2)
%         block_order{2} = '3'; % medical first
%     else
%         block_order{2} = '4'; % monetary first
%     end
% end
% 
% monstart = Datamon.trialTime(2*tpb+1).trialStartTime(4:6);
% medstart = Datamed.trialTime(2*tpb+1).trialStartTime(4:6);
% 
% if is_mon % domain is monetary
%     if monstart(1)*60 + monstart(2) > medstart(1)*60 + medstart(2) %time in minutes
%         block_order{3} = '6'; % medical first
%     else
%         block_order{3} = '5'; % monetary first
%     end
% else % domain is medical
%     if monstart(1)*60 + monstart(2) > medstart(1)*60 + medstart(2)
%         block_order{3} = '5'; % medical first
%     else
%         block_order{3} = '6'; % monetary first
%     end
% end
% 
% 
% monstart = Datamon.trialTime(3*tpb+1).trialStartTime(4:6);
% medstart = Datamed.trialTime(3*tpb+1).trialStartTime(4:6);
% 
% if is_mon % domain is monetary
%     if monstart(1)*60 + monstart(2) > medstart(1)*60 + medstart(2) %time in minutes
%         block_order{4} = '8'; % medical first
%     else
%         block_order{4} = '7'; % monetary first
%     end
% else % domain is medical
%     if monstart(1)*60 + monstart(2) > medstart(1)*60 + medstart(2)
%         block_order{4} = '7'; % medical first
%     else
%         block_order{4} = '8'; % monetary first
%     end
% end

end
