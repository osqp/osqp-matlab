function osqp_block_init(parentBlock)

    %------------------------------------------------------
    % Extract DialogParameters from the block
    refactor_style              = get_param(parentBlock,'refactor_style');
    warmstart_style             = get_param(parentBlock,'warmstart_style');
    has_external_q_input        = get_param(parentBlock,'has_external_q_input');
    has_external_bound_input    = get_param(parentBlock,'has_external_bound_input');
    
       
    %I want the ports to always appear in the same order
    %on the block if the same dialog values are selected.
    %this only seems possible if I manually track the
    %and assign the port numbers 
    portcount = 0;   %no ports assigned yet


    %------------------------------------------------------
    %Handle the configurable input ports that use pulldowns
    %------------------------------------------------------

    thePortnames         = {'warmstart_ext','refactor_ext'};
    thePulldownValues    = {warmstart_style, refactor_style};
 

    for i = 1:length(thePortnames)

        portname = thePortnames{i};
        pdVal    = thePulldownValues{i};

        %Configure the block with name = <portname> to be
        %of either constant or inport type
        %------------------------------------
        switch pdVal

            case {'Never','Always'}
                portcount = osqp_replace_block(parentBlock,portname,'Constant',portcount);

            case {'Triggered'}        %From external signal
                portcount = osqp_replace_block(parentBlock,portname,'Inport',portcount);

            otherwise
                error('Unrecognized pulldown value');

        end % switch

        %For constant cases, set an appropriate disable/enable signal
        %------------------------------------
        switch pdVal

            case 'Never'
                set_param([parentBlock, '/' portname],'Value','0')

            case 'Always'
                set_param([parentBlock, '/' portname],'Value','1')

            case 'Triggered'
                %do nothing, value is from external port

            otherwise
                error('Unrecognized pulldown value');

        end %switch

    end %end loop

    %------------------------------------------------------
    %Enable/Disable external matrix data if refactors are enabled
    %------------------------------------------------------

    thePortnames         = {'P_idx_ext','P_val_ext','A_idx_ext','A_val_ext'};

    %Configure the block with name = <portname> to be
    %of either constant or inport type
    %------------------------------------
    for j = 1:length(thePortnames)
        portname = thePortnames{j};

        switch refactor_style

            case {'Always','Triggered'}
                portcount = osqp_replace_block(parentBlock,portname,'Inport',portcount);

            case 'Never'
                portcount = osqp_replace_block(parentBlock,portname,'Constant',portcount);
                set_param([parentBlock, '/' portname],'Value','NaN')

            otherwise
                error('Unrecognized pulldown value');

        end


    end %end loop


    %------------------------------------------------------
    %Handle the configurable input ports that use tickboxes
    %------------------------------------------------------

    thePortGroups    = {{'q_ext'},{'l_ext','u_ext'}};
    theTickValues    = {has_external_q_input,has_external_bound_input};

    for i = 1:length(theTickValues)

        portnames = thePortGroups{i};
        tickVal    = theTickValues{i};

        %Configure the block with name = <portname> to be
        %of either constant or inport type
        %------------------------------------
        for j = 1:length(portnames)
            portname = portnames{j};

            switch tickVal

                case 'off'   %unticked
                    portcount = osqp_replace_block(parentBlock,portname,'Constant',portcount);
                    set_param([parentBlock, '/' portname],'Value','NaN')

                case 'on'   %ticked
                    portcount = osqp_replace_block(parentBlock,portname,'Inport',portcount);

                otherwise
                    error('Unrecognized dialog value');

            end %switch

        end % for j

    end %end loop

end %end function








function portcount = osqp_replace_block(parentBlock,portname,newtype,portcount)

    %replace the named block with one of the specified type, but
    %only do this if the new type is different from the current one

    %this is necessary since replacing an input port block with
    %the same type has the effect of disconnecting the port
    
    blockname = [parentBlock, '/' portname];
    
    oldtype = get_param(blockname,'BlockType');

    if(~strcmp(oldtype,newtype))
        replace_block(blockname,'Name',portname,newtype,'noprompt');
    end
    
    %if this is an Inport, assign it the next port number and increment
    if(strcmp(newtype,'Inport'))
        portcount = portcount + 1;
        set_param(blockname,'Port',num2str(portcount));
    end
        
    
   

end %end function
