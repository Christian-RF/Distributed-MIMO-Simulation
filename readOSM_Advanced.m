function buildings = readOSM_Advanced(filename)
    % READOSM_ADVANCED Reads OSM XML, extracting Lat, Lon, and Height of buildings.
    fprintf('Parsing OSM file: %s ...\n', filename);
    try
        xDoc = xmlread(filename);
    catch
        error('Could not read file.');
    end
    
    % Index Nodes 
    allNodes = xDoc.getElementsByTagName('node');
    nodeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 0:allNodes.getLength-1
        item = allNodes.item(i);
        id = char(item.getAttribute('id'));
        coords = [str2double(item.getAttribute('lat')), str2double(item.getAttribute('lon'))];
        nodeMap(id) = coords;
    end
    
    % Index Ways 
    allWays = xDoc.getElementsByTagName('way');
    wayMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    fprintf('Indexing %d ways...\n', allWays.getLength);
    
    for i = 0:allWays.getLength-1
        thisWay = allWays.item(i);
        wayID = char(thisWay.getAttribute('id'));
        
        ndRefs = thisWay.getElementsByTagName('nd');
        wLat = []; wLon = [];
        
        for k = 0:ndRefs.getLength-1
            ref = char(ndRefs.item(k).getAttribute('ref'));
            if isKey(nodeMap, ref)
                c = nodeMap(ref);
                wLat(end+1) = c(1); %#ok<AGROW>
                wLon(end+1) = c(2); %#ok<AGROW>
            end
        end
        
        if ~isempty(wLat)
            s.Lat = wLat; s.Lon = wLon;
            s.Tags = getTags(thisWay); 
            wayMap(wayID) = s;
        end
    end
    
    % Extract Buildings 
    buildings = struct('Lat', {}, 'Lon', {}, 'Height', {}); % Added Height field
    bCount = 0;
    
    % Process Simple Ways
    wayKeys = keys(wayMap);
    for i = 1:length(wayKeys)
        w = wayMap(wayKeys{i});
        if isKey(w.Tags, 'building')
            bCount = bCount + 1;
            buildings(bCount).Lat = w.Lat;
            buildings(bCount).Lon = w.Lon;
            % Extract Height
            buildings(bCount).Height = parseHeight(w.Tags);
        end
    end
    
    % Process Relations
    allRelations = xDoc.getElementsByTagName('relation');
    fprintf('Checking %d relations...\n', allRelations.getLength);
    
    for i = 0:allRelations.getLength-1
        thisRel = allRelations.item(i);
        tags = getTags(thisRel);
        
        if isKey(tags, 'building') || (isKey(tags, 'type') && strcmp(tags('type'), 'multipolygon') && isKey(tags, 'building'))
             % Use relation tags for height if available, otherwise default
            relHeight = parseHeight(tags);
            
            members = thisRel.getElementsByTagName('member');
            for k = 0:members.getLength-1
                mem = members.item(k);
                type = char(mem.getAttribute('type'));
                ref = char(mem.getAttribute('ref'));
                role = char(mem.getAttribute('role'));
                
                if strcmp(type, 'way') && (strcmp(role, 'outer') || strcmp(role, ''))
                    if isKey(wayMap, ref)
                        part = wayMap(ref);
                        bCount = bCount + 1;
                        buildings(bCount).Lat = part.Lat;
                        buildings(bCount).Lon = part.Lon;
                        
                        % If the way itself has height, use it, otherwise use Relation height
                        partHeight = parseHeight(part.Tags);
                        if partHeight > 10 % If way has specific height
                            buildings(bCount).Height = partHeight;
                        else
                            buildings(bCount).Height = relHeight;
                        end
                    end
                end
            end
        end
    end
    fprintf('Found %d building parts.\n', bCount);
end

% Helpers 
function tMap = getTags(elem)
    tMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
    tags = elem.getElementsByTagName('tag');
    for k = 0:tags.getLength-1
        t = tags.item(k);
        key = char(t.getAttribute('k'));
        val = char(t.getAttribute('v'));
        tMap(key) = val;
    end
end

function h = parseHeight(tags)
    % Logic to find height or levels, or return default
    if isKey(tags, 'height')
        % Sometimes height is "12 m", strip unit
        raw = tags('height');
        val = str2double(regexprep(raw, '[^0-9.]', '')); 
        if ~isnan(val), h = val; return; end
    end
    
    if isKey(tags, 'building:levels')
        val = str2double(tags('building:levels'));
        if ~isnan(val)
            h = val * 3.5; % Assume 3.5m per floor
            return; 
        end
    end
    
    h = 10; % Default fallback height (meters)
end