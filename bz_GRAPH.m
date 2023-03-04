classdef bz_GRAPH
    %BZ_GRAPH band structure for Orthorhombic crystal
    %   Detailed explanation goes here
    
    properties
        ICSD;
        TQC_profile;
        bandStructure;
        energyDensityStructure;

        TQC_profile_nosoc;
        bandStructure_nosoc;
        energyDensityStructure_nosoc;

        kPath;

        kPathGraph;
        kPathGraph_nosoc;

    end
    
    methods
        
        function obj = bz_GRAPH(ICSD)
            %BZ_GRAPH Construct an instance of this class
            %   Detailed explanation goes here
            obj.ICSD = ICSD;
            url = sprintf('https://www.topologicalquantumchemistry.com/api/v1/compounds/icsd=%d/',ICSD);
            obj.TQC_profile = webread(url);

            url2 = sprintf('https://www.topologicalquantumchemistry.com/api/v1/compounds/%d/bandStructure/',obj.TQC_profile.id);
            obj.bandStructure = webread(url2);

            url3 = sprintf('https://www.topologicalquantumchemistry.com/api/v1/compounds/%d/energyDensityStructure/',obj.TQC_profile.id);
            obj.energyDensityStructure = webread(url3);

            url = sprintf('https://www.topologicalquantumchemistry.com/api/v1/compounds/icsd=%d/nosoc/',ICSD);
            obj.TQC_profile_nosoc = webread(url);

            url2 = sprintf('https://www.topologicalquantumchemistry.com/api/v1/compounds/%d/nosoc/bandStructure/',obj.TQC_profile_nosoc.id);
            obj.bandStructure_nosoc = webread(url2);

            url3 = sprintf('https://www.topologicalquantumchemistry.com/api/v1/compounds/%d/nosoc/energyDensityStructure/',obj.TQC_profile_nosoc.id);
            obj.energyDensityStructure_nosoc = webread(url3);


            obj.kPath = obj.getKPath();

            obj.kPathGraph = obj.praseBandStructure2(obj.bandStructure);
            obj.kPathGraph_nosoc = obj.praseBandStructure2(obj.bandStructure_nosoc);


        end

        function kPath = getKPath(obj)
            % get k-path info from SeeK-Path website

            % save .cif file
            if ~exist('./cifFile', 'dir')
                mkdir('./cifFile')
            end
            websave(['./cifFile/' num2str(obj.ICSD) '.cif'], obj.TQC_profile.cifFile);

            % run python code
            pyrunfile(append("cif2json.py './cifFile/", num2str(obj.ICSD), ".cif'"))

            % prase json
            fname = ['./cifFile/', num2str(obj.ICSD), '.json']; 
            fid = fopen(fname); 
            raw = fread(fid,inf); 
            str = char(raw'); 
            fclose(fid); 
            kPath = jsondecode(str);

        end

        function KPG = praseBandStructure2(obj,bs)
            
            % length of total k
            l = length(bs.kp);
            % total k
            k = zeros(1,l);
            for i = 1:l
                k(i) = bs.kp(i).v;
            end
            
            % number of bands
            n = length(bs.eb);
            % energy of all bands
            eb = zeros(n,l);
            for m = 1:l
                for n = 1:n
                    eb(n,m) = bs.eb(n).i(m).e;
                end
            end

            % construct the graph
            startFlag = 1;
            name = {};
            s = [];
            t = [];
            s_ind = [];
            t_ind = [];
            
            for i = 1:l % for every k
                if bs.kp(i).h == 1 % if it has high symmetry
            
                    if startFlag == 1 % if it is the beginning
            
                        hn = bs.kp(i).hn; % high symmetry name
                        if contains(bs.kp(i).hn,'|')
                            temp = split(bs.kp(i).hn,'|');
                            hn = temp{2};
                        end
                
                        % add high symmetry point to the name list
                        if ~any(strcmp(name,hn))
                            name{end+1} = hn;
                        end
                        nodeIndex = find(strcmp(name,hn));
            
                        s(end+1) = nodeIndex;
                        s_ind(end+1) = i;
                        startFlag = 0;

                    else % if it is the ending
            
                        hn = bs.kp(i).hn;
                        if contains(bs.kp(i).hn,'|')
                            temp = split(bs.kp(i).hn,'|');
                            hn = temp{1};
                        end
            
                        if ~any(strcmp(name,hn))
                            name{end+1} = hn;
                        end
                        nodeIndex = find(strcmp(name,hn));
            
                        t(end+1) = nodeIndex;
                        t_ind(end+1) = i;
                        startFlag = 1;
                    end
                end
            end
            
            % results
            edge_num = length(s);
            k_list = cell(edge_num,1);
            kl_list = cell(edge_num,1);
            eb_list = cell(edge_num,1);

            range_list = cell(edge_num,2);

            for i = 1:edge_num
                k_list{i} = 2*pi*(k(s_ind(i):t_ind(i)) - k(s_ind(i)));
                kl_list{i} = num2str(2*pi*(k(t_ind(i)) - k(s_ind(i))));
                eb_i = eb(:,s_ind(i):t_ind(i));
                if s(i)>t(i)
                    eb_i = flip(eb_i,2);
                end
                eb_list{i} = eb_i;
                range_list{i,1} = min(eb_i,[],2);
                range_list{i,2} = max(eb_i,[],2);
            end
            
            % specify position of each high symmetry point
            node_num = length(name);
            hp_position = cell(node_num,1);
            for i=1:node_num
                hn = name{i};
                if hn == 'G'
                    hn = 'GAMMA';
                end
                hp_position{i} = obj.kPath.kpoints.(hn);
            end
            


            EdgeTable = table([s' t'],k_list,kl_list,eb_list,range_list, 'VariableNames',{'EndNodes' 'K_axis' 'K_length' 'BandEnergy' 'Range'});
%             NodeTable = table(name',hp_position,'VariableNames',{'Name' 'Position'});
            NodeTable = table(name',hp_position,'VariableNames',{'Name','Position'});
            KPG = graph(EdgeTable,NodeTable);

        end

        % ---- plot -------------------------------
        function plotUnitCell(obj)
            a1 = obj.kPath.primitive_lattice(1,:);
            a2 = obj.kPath.primitive_lattice(2,:);
            a3 = obj.kPath.primitive_lattice(3,:);

            figure;
            hold on
            plot3([0 a1(1)],[0 a1(2)],[0 a1(3)]);
            text(a1(1),a1(2),a1(3),"a_1");
            plot3([0 a2(1)],[0 a2(2)],[0 a2(3)]);
            text(a2(1),a2(2),a2(3),"a_2");
            plot3([0 a3(1)],[0 a3(2)],[0 a3(3)]);
            text(a3(1),a3(2),a3(3),"a_3");

            x = [a1(1) -a1(1) -a1(1) a1(1) a1(1); a1(1) -a1(1) -a1(1) a1(1) a1(1)]/2;
            y = [a2(2) a2(2) -a2(2) -a2(2) a2(2); a2(2) a2(2) -a2(2) -a2(2) a2(2)]/2;
            z = [a3(3)*ones(1,5); -a3(3)*ones(1,5)]/2;
            surf(x, y, z,'FaceColor','w','FaceAlpha',0.1);
            patch(x', y', z', 'b','FaceAlpha',0.1);
            set(gca,'DataAspectRatio',[1 1 1]);

            for i=1:length(obj.kPath.primitive_symbols)
                x = obj.kPath.primitive_positions_cartesian(i,1);
                y = obj.kPath.primitive_positions_cartesian(i,2);
                z = obj.kPath.primitive_positions_cartesian(i,3);
                symb = obj.kPath.primitive_symbols{i};
                scatter3(x,y,z);
                text(x,y,z,symb);
            end

        end

        function plotBrillouinZone(obj)
            a1 = obj.kPath.b1;
            a2 = obj.kPath.b2;
            a3 = obj.kPath.b3;

            figure;
            hold on
            plot3([0 a1(1)],[0 a1(2)],[0 a1(3)]);
            text(a1(1),a1(2),a1(3),"b_1");
            plot3([0 a2(1)],[0 a2(2)],[0 a2(3)]);
            text(a2(1),a2(2),a2(3),"b_2");
            plot3([0 a3(1)],[0 a3(2)],[0 a3(3)]);
            text(a3(1),a3(2),a3(3),"b_3");

            x = [a1(1) -a1(1) -a1(1) a1(1) a1(1); a1(1) -a1(1) -a1(1) a1(1) a1(1)]/2;
            y = [a2(2) a2(2) -a2(2) -a2(2) a2(2); a2(2) a2(2) -a2(2) -a2(2) a2(2)]/2;
            z = [a3(3)*ones(1,5); -a3(3)*ones(1,5)]/2;
            surf(x, y, z,'FaceColor','w','FaceAlpha',0.1);
            patch(x', y', z', 'b','FaceAlpha',0.1);
            set(gca,'DataAspectRatio',[1 1 1]);

            fn = fieldnames(obj.kPath.kpoints);
            for i=1:length(fn)
                POS = obj.kPath.kpoints.(fn{i});
                x = POS(1);
                y = POS(2);
                z = POS(3);

                symb = fn{i};
                scatter3(x,y,z);
                text(x,y,z,symb);
            end

            for i=1:length(obj.kPath.path)
                S = obj.kPath.path{i};
                X1 = obj.kPath.kpoints.(S{1});
                X2 = obj.kPath.kpoints.(S{2});
                plot3([X1(1) X2(1)],[X1(2) X2(2)],[X1(3) X2(3)]);
            end
        end

        function plotGraph(obj)
            figure
            plot(obj.kPathGraph,'NodeLabel',obj.kPathGraph.Nodes.Name,'EdgeLabel',obj.kPathGraph.Edges.K_length);

        end
        
        function plotBZonCurrentAxes(obj,vl,hl)

            axes = gca;
            % get axes limits
            x_lim = axes.XLim;
            y_lim = axes.YLim;

            switch vl
                case 'a'
                    bx = obj.kPath.b1(1);
                case 'b'
                    bx = obj.kPath.b2(2);
                case 'c'
                    bx = obj.kPath.b3(3);
                otherwise
                    bx = 0;
            end

            switch hl
                case 'a'
                    by = obj.kPath.b1(1);
                case 'b'
                    by = obj.kPath.b2(2);
                case 'c'
                    by = obj.kPath.b3(3);
                otherwise
                    by = 0;
            end

            x_n2(1) = round(x_lim(1)/bx-0.5)+0.5;
            x_n2(2) = round(x_lim(2)/bx-0.5)+0.5;
            y_n2(1) = round(y_lim(1)/by-0.5)+0.5;
            y_n2(2) = round(y_lim(2)/by-0.5)+0.5;

            % plot energy contour lines
            lines = {};
            hold(axes,'on');
            if bx~= 0
                for xx = x_n2(1):x_n2(2)
                    lines{end+1} = xline(xx*bx,'--','Color',[0.6350 0.0780 0.1840],'Alpha',0.5,'LineWidth',1);
                end
            end

            if by~= 0
                for yy = y_n2(1):y_n2(2)
                    lines{end+1} = yline(yy*by,'--','Color',[0.6350 0.0780 0.1840],'Alpha',0.5,'LineWidth',1);
                end
            end
            hold(axes,'off');
        end

        % ---- band -----------------------------
        function [K,EB] = getKnE(obj,hns,hnt,soc)
            if strcmp(soc,'soc')
                GR = obj.kPathGraph;
            elseif strcmp(soc,'nosoc')
                GR = obj.kPathGraph_nosoc;
            else
                disp('soc/nosoc?');
                return
            end
            % get K in 1d and energy of bands
            edgeIndex = findedge(GR,hns,hnt);
            if edgeIndex == 0
                disp('Not found.');
                return
            end
            K = GR.Edges.K_axis{edgeIndex};
            EB = GR.Edges.BandEnergy{edgeIndex};
            
            % flip the direction if the direction is inversed
            if strcmp(hns,GR.Edges.EndNodes{edgeIndex,2})
                EB = flip(EB,2);
            end
        end
        
        function plotBand(obj,hns,hnt,soc)

            [k,E] = obj.getKnE(hns,hnt,soc);

            figure
            plot(k,E,'Color',"#A2142F",'LineWidth',1);
            title(append('Dispersion along ',hns,'-',hnt));
            xlim([min(k) max(k)]);
            ylim([-5 1]);
            xlabel('K [Å^{-1}]');
            ylabel('Binding Energy [eV]');
        end

        function plotBandOnCurrentAxes(obj,hns,hnt,soc)
            
            [k,E] = obj.getKnE(hns,hnt,soc);

            if strcmp(soc,'soc')
                clr = [0.8500 0.3250 0.0980 0.3];
            elseif strcmp(soc,'nosoc')
                clr = [0.4940 0.1840 0.5560 0.3];
            else
                disp('soc/nosoc?');
                return
            end

            % get current axes
            axes = gca;

            % get axes limits
            LL = axes.XLim;
            KM = max(k);

            % plot bands
            hold(axes,'on');
            for i = round(LL(1)/KM/2):round(LL(2)/KM/2)
                plot(axes,k+i*KM*2,E,'Color',clr,'LineWidth',1);
                plot(axes,-k+i*KM*2,E,'Color',clr,'LineWidth',1);
            end
            % plot fermi level
            yline(axes,0,'--b','LineWidth',1,'Alpha',0.5);
            hold(axes,'off');
        end

        % ---- get --------------------------------
        function [K3,EB] = getK3nE(obj,hns,hnt,soc)
            if strcmp(soc,'soc')
                GR = obj.kPathGraph;
            elseif strcmp(soc,'nosoc')
                GR = obj.kPathGraph_nosoc;
            else
                disp('soc/nosoc?');
                return
            end

            % get K in 3d space and energy of bands

            % get energy of bands
            edgeIndex = findedge(GR,hns,hnt);
            if edgeIndex == 0
                disp('Not found.');
                return
            end
            K = GR.Edges.K_axis{edgeIndex};
            EB = GR.Edges.BandEnergy{edgeIndex};
            % flip the direction if the direction is inversed
            if strcmp(hns,GR.Edges.EndNodes{edgeIndex,2})
                EB = flip(EB,2);
            end

            % get K in 3d space
            S_nodeIndex = findnode(GR,hns);
            S_position = GR.Nodes.Position{S_nodeIndex};
            T_nodeIndex = findnode(GR,hnt);
            T_position = GR.Nodes.Position{T_nodeIndex};

            K3(1,:) = linspace(S_position(1),T_position(1),length(K));
            K3(2,:) = linspace(S_position(2),T_position(2),length(K));
            K3(3,:) = linspace(S_position(3),T_position(3),length(K));

        end

        function plot_BZ_test(obj,hn_list,Emin,Emax)

            K3_l = [];
            EB_l = [];

            nbr_hn = length(hn_list);
            nodeIndex = findnode(obj.kPathGraph,hn_list);

            for i=1:nbr_hn-1
                hns = hn_list{i};
                hnt = hn_list{i+1};

                [K3_i,EB_i] = obj.getK3nE(hns,hnt);

                K3_l = horzcat(K3_l,K3_i);
                EB_l = horzcat(EB_l,EB_i);

            end

            % find n by energy range
            n = [];
            for i = 1:size(EB_l,1)
                mi = min(EB_l(i,:));
                ma = max(EB_l(i,:));
                if mi < Emax && ma > Emin
                    n(end+1) = i;
                end
            end


            figure;
            hold on
            plot3(K3_l(1,:), K3_l(2,:),EB_l(n,:),'LineWidth',1.5);
            plot3(K3_l(2,:), K3_l(1,:),EB_l(n,:),'LineWidth',1.5);
            plot3(-K3_l(1,:), K3_l(2,:),EB_l(n,:),'LineWidth',1.5);
            plot3(K3_l(2,:), -K3_l(1,:),EB_l(n,:),'LineWidth',1.5);
            plot3(K3_l(1,:), -K3_l(2,:)+obj.kPath.b1(1),EB_l(n,:),'LineWidth',1.5);
            plot3(-K3_l(2,:)+obj.kPath.b1(1),K3_l(1,:),EB_l(n,:),'LineWidth',1.5);
%             plot3(K3_l(2,:), -K3_l(1,:)+obj.kPath.b1(1),EB_l(n,:),'LineWidth',1.5);

        end

        function [K3, EB] = getK3nE_l(obj,hn_list,soc)
            if strcmp(soc,'soc')
                GR = obj.kPathGraph;
            elseif strcmp(soc,'nosoc')
                GR = obj.kPathGraph_nosoc;
            else
                disp('soc/nosoc?');
                return
            end
            
            K3 = [];
            EB = [];
            nbr_hn = length(hn_list);
            nodeIndex = findnode(GR,hn_list);

            for i=1:nbr_hn-1
                hns = hn_list{i};
                hnt = hn_list{i+1};

                edgeIndex = findedge(GR,hns,hnt);
                if edgeIndex == 0
                    disp('Not found.');
                    return
                end

                EB_i = GR.Edges.BandEnergy{edgeIndex};
                if strcmp(hns,GR.Edges.EndNodes{edgeIndex,2})
                    EB_i = flip(EB_i,2);
                end

                n = linspace(0,1,size(EB_i,2));

                posi_s = GR.Nodes.Position{nodeIndex(i)};
                posi_t = GR.Nodes.Position{nodeIndex(i+1)};
                K3_i = posi_s + (posi_t - posi_s)*n;

                K3 = horzcat(K3,K3_i);
                EB = horzcat(EB,EB_i);
            end
        end
    end
end

