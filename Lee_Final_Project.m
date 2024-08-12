function Lee_Final_Project
%   This code simulates the spread of COVID-19 within a 2D box. 
%   Agents move in the box at a fixed speed along their respective directions, 
%   slightly altering their direction at each step. For each timestep 
%   (considered an hour), agents wander the grid. When an agent enters a 
%   grid cell containing another agent, it faces a probability of infection 
%   if one of the agents is already infected. Once an agent has been 
%   infected for six days, their status changes to 'Quarantine,' 
%   and they remain stationary in the grid cell until their quarantine 
%   period ends. During quarantine, agents still have a 0.001 probability 
%   of spreading the virus. Upon completing the quarantine period, 
%   agents become immune.
%   

% Agent Parameters
N = 500; % initial number of agents
v = 0.0125;  % speed of agents (slow: 0.0125, normal: 0.025, fast: 0.05)

% Time Parameters
hours_in_day = 24; % consider each step an hour
ns = 115 * hours_in_day; % number of steps (hours)
max_infection = 6 * hours_in_day; % duration of infection before quarantine
quarantine_over = 20 * hours_in_day; % duration after which agent is no longer infectious, and is immune

% Domain Parameters
ax = [0 10 0 10]; % domain of the box
cx = ax(1):ax(2); % cell boundary(x)
cy = ax(3):ax(4); % cell boundary(y)

% Initialization
X = [ax(1)+(ax(2)-ax(1))*rand(N,1), ax(3)+(ax(4)-ax(3))*rand(N,1)]; % Initial positions of agents (spread throughout the entire domain)
I = double(rand(N,1)<0.019); % infection status for each agent (0=negative, >0 is the number of hours that they have been positive)
I = I .* floor(rand(N,1) * (max_infection / hours_in_day)); % give each agent that is initially positive a random duration of infection (under the amount before they quarantine)
D = rand(N,1)*2*pi; % initial angles of direction of agent
H = zeros(0, 5); % array to hold the history of agents (time, infected, immune, quarantined, susceptible)
M = zeros(N,1); % array to hold the immune status for each agent (0 = not immune, 1 = immune)

% Calculate the baseline probability for an agent sharing a square with an
% infected agent to become infected
r0 = 3.0; % amount of agents a single infected agent will infect(in range 3.28 and 2.79)
avg_pop_density = N / ((ax(2)-ax(1)) * (ax(4)-ax(3))); % average amount of people in cell to infect
avg_cells_visited = (max_infection * v) / 4; % Average number of cells to visit, divided by 4
p_base = r0 / (avg_pop_density * avg_cells_visited * max_infection); % max_infection is amount of opportunities to infect

for j = 1:ns % loop over steps
    % Update positions and directions
    for k = 1:N
        if I(k) >= max_infection && I(k) < quarantine_over % if the agent is in quarantine
            continue; % Skip updating the position

        elseif I(k) < max_infection % if not in quarantine,
            X(k,1) = X(k,1) + v*cos(D(k)); % Update x-position
            X(k,2) = X(k,2) + v*sin(D(k)); % Update y-position
            D(k) = D(k) + 0.1*randn; % change direction of motion
        end
    end
    
    % Contain agents in the domain: let agents "bounce" off walls
    ind = (X(:,1)<ax(1)&cos(D)<0)|... % Agents hitting a boundary
        (X(:,1)>ax(2)&cos(D)>0); % horizontally
    D(ind) = pi-D(ind); % reverse x-direction
    ind = (X(:,2)<ax(3)&sin(D)<0)|... % Agents hitting a boundary
        (X(:,2)>ax(4)&sin(D)>0); % vertically
    D(ind) = -D(ind); % reverse y-direction
    X(:,1) = min(max(X(:,1),ax(1)),ax(2)); % move agents outside of
    X(:,2) = min(max(X(:,2),ax(3)),ax(4)); % domain back onto boundary

    [occupied_cells,~,agents_in_cell] = unique(floor(X),'rows'); 
    cells_with_multiple_agents = occupied_cells(accumarray(agents_in_cell,1)>=2,:); % define cells with multiple agents

    for k = 1:size(cells_with_multiple_agents,1) % for every cell with adequate number of agents
        agents_in_this_cell = zeros(0,0);
        quarantine_present = 0; % Number of the agents in quarantine in the cell
        covid_present = 0; % Number of the agents who are infected in the cell

        for l = 1:size(X,1) % Check if they are in the same cell add them to a list and check to see if they have COVID
            if (X(l,1) > cells_with_multiple_agents(k,1) && X(l,1) - cells_with_multiple_agents(k,1) <= 1 && X(l,2) > cells_with_multiple_agents(k,2) && X(l,2) - cells_with_multiple_agents(k,2) <= 1) % to see if they are in a cell together
                agents_in_this_cell = [agents_in_this_cell;l];
                
                if I(l) >= 1 && I(l) <= max_infection && M(l) ~= 1
                    covid_present = 1; % Check if any agents have COVID (and have yet to quarantine)
                
                elseif I(l) > max_infection && M(l) ~= 1 % if the agent quarantine status and not immune
                    quarantine_present = 1; % Increase the number of quarantine agent
                end
            end
        end

        if covid_present % If COVID is present in this cell
            
            for l = 1:size(agents_in_this_cell) % find the agents without COVID and with a probability of the base COVID probability times the adjusted odds ratio 
                if I(agents_in_this_cell(l)) == 0 && rand < p_base
                    I(agents_in_this_cell(l)) = 1; % give them COVID (and start the countdown before their quarantine)
                end
            end
        elseif quarantine_present % If quarantine agent is present in this cell
            for l = 1:size(agents_in_this_cell) 
                if I(agents_in_this_cell(l)) == 0 && rand < 0.001
                    I(agents_in_this_cell(l)) = 1; 
                end
            end
        end
    end

    for k = 1:size(I) % check every agent
        if I(k) >= quarantine_over % if quarantine is over 
            M(k) = 1; % set agent is immune
            I(k) = 0; % set the agent recovered
        elseif I(k)>=1 % to see if they are infected
            I(k) = I(k) + 1; % if so, add a day to their time before quarantine
        end
    end

    current_infected = nnz(I > 0 & I < max_infection); % count the number of infected agents
    current_immune = nnz(I == 0 & M == 1);  % count the number of immune agents
    current_quarantined = nnz(I >= max_infection & M == 0);  % count the number of quarantined agents
    current_susceptible = nnz(I==0 & M==0);  % count the number of susceptible agents
    H(end+1, :) = [j, current_infected, current_immune, current_quarantined, current_susceptible]; % Add the number of infected, immuned, quaratined, and susceptible agents
    
    % Plotting
    if floor(j/12)*12==j % plot for every 12 hours

        % Subplot 1:
        % clf
        % subplot(1,2,1);
        % plot([1;1]*cx,ax([3,4])'*(cx*0+1),'k-',... % draw boundaries
        %     ax([1,2])'*(cy*0+1),[1;1]*cy,'k-') % of cells
        % hold on
        % 
        % leg = legend("");
        % 
        % leg.Location = 'northeastoutside'; % set location of legend 
        % leg.FontSize = 12; % Adjust the font size as needed
        % 
        % 
        % % Plot points (agents)
        % 
        % % Susceptible agents
        % ind = I==0 & M==0;
        % plot(X(ind,1),X(ind,2), '.','markersize',12,'color',[0 0 0.5], 'DisplayName', sprintf('Susceptible (%d)',nnz(ind)));    
        % % infected agents 
        % ind = I >0 & I < max_infection; % infected agents that have not yet quarantined
        % plot(X(ind,1),X(ind,2),'.','markersize',12,'color',[0.9 0 0.1], 'DisplayName', sprintf('Infected (%d)',nnz(ind)))
        % % quarantined
        % ind = I >= max_infection & M == 0; 
        % plot(X(ind,1),X(ind,2),'o','markersize',6,'color',[0.8500 0.3250 0.0980], 'DisplayName', sprintf('Quarantined (%d)',nnz(ind)))
        % % Plot immune agents
        % ind = I == 0 & M==1;
        % plot(X(ind,1),X(ind,2),'.','markersize',12,'color',[0 0.5 0], 'DisplayName', sprintf('Immune (%d)',nnz(ind)))
        % 
        % hold off
        % axis equal xy, axis(ax)
        % title('Spread of COVID')
        % subtitle(sprintf('(%d agents, %d days and %d hours)', current_susceptible, (floor(j/24)), mod(j,24)))
        % 
        % % Subplot 2
        % subplot(1,2,2);
        % hold on
        % plot(H(:,1), H(:,5), '-', 'color',[0 0 0.5], 'LineWidth', 2); % graph to show how many agents are susceptible
        % plot(H(:,1), H(:,2), '-', 'color',[0.9 0 0.1], 'LineWidth', 2); % graph to show how many agents are positive for COVID
        % plot(H(:,1), H(:,4), '-', 'color',[0.8500 0.3250 0.0980], 'LineWidth', 2); % graph to show how many agents are quarantined
        % plot(H(:,1), H(:,3), '-', 'color',[0 0.5 0], 'LineWidth', 2); % graph to show how many agents are immune
        % text(1,current_infected, sprintf("Current infected agent size (%d)", current_infected),'color', '#555555') 
        % yline(current_infected, ':','color', '#555555'); % line to show how many agents are not quarantining
        % axis([1 max(j,10) 0 max(max(H(:,2:end)))]); % Set y-axis limit based on the maximum value across all columns of H except the first one
        % xlabel('Hours'), ylabel('Number of Agents')
        % legend('Susceptible','Infected','Quarantined', 'Immune','FontSize',12)
        % 
        % drawnow

        if current_infected == 0 && current_quarantined == 0 % If there are no more infected and quarantined agents, exit the loop
            break; % Exit the loop
        end

    end
end

% Final graph
figure
hold on
% plot(H(:,1), H(:,5), '-', 'color',[0 0 0.5],'LineWidth', 2); % graph to show how many agents are susceptible
plot(H(:,1), H(:,2), '-', 'color',[0.9 0 0.1],'LineWidth', 2); % graph to show how many agents are positive for COVID
% plot(H(:,1), H(:,4), '-', 'color',[0.8500 0.3250 0.0980],'LineWidth', 2); % graph to show how many agents are quarantined
% plot(H(:,1), H(:,3), '-', 'color',[0 0.5 0],'LineWidth', 2); % graph to show how many agents are immune

% Find the index of the maximum number of infected agents
[~, max_infected_idx] = max(H(:,2));

% Plot the maximum infection marker
plot(H(max_infected_idx, 1), H(max_infected_idx, 2), 'ro', 'MarkerSize', 12, 'DisplayName', 'Max Infection');

% Add text label showing the number of infected agents
text(H(max_infected_idx, 1), H(max_infected_idx, 2), sprintf('Max Infected: %d', H(max_infected_idx, 2)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left','FontSize',12);
axis([1 max(j,10) 0 max(max(H(:,2:end)))]); % Set y-axis limit based on the maximum value across all columns of H except the first one

xlabel('Hours'), ylabel('Number of Agents')
legend('Susceptible','Infected','Quarantined', 'Immune')
saveas(gcf, 'graph.png');
