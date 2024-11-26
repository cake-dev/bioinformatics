// Default dimensions - increased for better visualization
const DEFAULT_WIDTH = 1200;
const DEFAULT_HEIGHT = 800;

document.getElementById("graphForm").addEventListener("submit", function (event) {
    event.preventDefault();

    const k = parseInt(document.getElementById("kValue").value);
    const dna = document.getElementById("dnaString").value.trim();
    if (k > dna.length) {
        alert("k must be smaller than or equal to the length of the DNA string.");
        return;
    }

    const { nodes, links, kmers } = buildDeBruijnGraph(dna, k);
    renderGraph(nodes, links, k, dna);
    renderKmerTable(kmers);
});

// BuildDeBruijnGraph function remains the same
function buildDeBruijnGraph(sequence, k) {
    const kmers = [];
    const links = [];
    const nodes = {};
    
    const startPad = '$'.repeat(k-1);
    const paddedSequence = startPad + sequence + '$';
    
    for (let i = 0; i <= paddedSequence.length - k; i++) {
        const kmer = paddedSequence.substring(i, i + k);
        const prefix = kmer.substring(0, k-1);
        const suffix = kmer.substring(1, k);
        const transition = kmer[k-1];
        
        kmers.push({
            number: kmers.length + 1,
            length: prefix === startPad ? 0 : 1,
            vertex: prefix,
            word: transition
        });

        links.push({
            source: prefix,
            target: suffix,
            transition: transition,
            originalIndex: i  // Store the original position in the sequence
        });
    }

    links.forEach(function(link) {
        link.source = nodes[link.source] || (nodes[link.source] = {name: link.source});
        link.target = nodes[link.target] || (nodes[link.target] = {name: link.target});
    });

    return { nodes: d3.values(nodes), links, kmers };
}

function findPathToMatch(nodes, links, originalString, k) {
    const adj = {};
    nodes.forEach(node => {
        adj[node.name] = [];
    });

    // Sort links by their original index to follow the sequence order
    links.forEach(link => {
        adj[link.source.name].push({
            target: link.target.name,
            transition: link.transition,
            originalIndex: link.originalIndex,
            used: false
        });
    });

    for (let node in adj) {
        adj[node].sort((a, b) => a.originalIndex - b.originalIndex);
    }

    const startNode = originalString.substring(0, k - 1); // Start at the first k-1 characters
    const path = [];
    let currentNode = startNode;

    while (true) {
        const edges = adj[currentNode];
        if (!edges || edges.length === 0) break;

        // Find the next unused edge
        const nextEdge = edges.find(edge => !edge.used);
        if (!nextEdge) break;

        nextEdge.used = true;
        path.push({
            from: currentNode,
            to: nextEdge.target,
            transition: nextEdge.transition
        });

        currentNode = nextEdge.target;
    }

    return path;
}


function findEulerianPath(nodes, links) {
    // Create adjacency list
    const adj = {};
    nodes.forEach(node => {
        adj[node.name] = [];
    });
    
    links.forEach(link => {
        adj[link.source.name].push({
            target: link.target.name,
            transition: link.transition
        });
    });

    // Find start node (node with one more outgoing edge)
    let startNode = nodes[0].name;
    for (let node in adj) {
        if (adj[node].length > 0) {
            startNode = node;
            break;
        }
    }

    // Perform Hierholzer's algorithm
    const path = [];
    const stack = [startNode];
    
    while (stack.length > 0) {
        const currentNode = stack[stack.length - 1];
        
        if (adj[currentNode].length > 0) {
            const nextEdge = adj[currentNode].pop();
            stack.push(nextEdge.target);
            path.push({
                from: currentNode,
                to: nextEdge.target,
                transition: nextEdge.transition
            });
        } else {
            stack.pop();
        }
    }

    return path;
}

function renderKmerTable(kmers) {
    const tableContainer = document.getElementById('kmer-table') || document.createElement('div');
    tableContainer.id = 'kmer-table';
    tableContainer.innerHTML = '';
    
    const table = document.createElement('table');
    
    // Create header
    const header = table.createTHead();
    const headerRow = header.insertRow();
    ['#', 'L', 'Ver', 'W'].forEach(text => {
        const th = document.createElement('th');
        th.textContent = text;
        headerRow.appendChild(th);
    });

    // Create body
    const tbody = table.createTBody();
    kmers.forEach((kmer, index) => {
        const row = tbody.insertRow();
        [
            index + 1,
            kmer.length,
            kmer.vertex,
            kmer.word
        ].forEach(text => {
            const td = document.createElement('td');
            td.textContent = text;
            row.appendChild(td);
        });
    });

    tableContainer.appendChild(table);
    document.querySelector('.container').appendChild(tableContainer);
}

function renderGraph(graphNodes, graphLinks, k, originalString) {
    // Create reconstruction display and controls
    const controlsContainer = document.createElement('div');
    controlsContainer.style.cssText = `
        padding: 20px;
        margin: 20px 0;
        display: flex;
        gap: 20px;
        align-items: center;
    `;
    
    const reconstructionDisplay = document.createElement('div');
    reconstructionDisplay.id = 'reconstruction-display';
    reconstructionDisplay.style.cssText = `
        flex-grow: 1;
        padding: 20px;
        font-family: monospace;
        font-size: 18px;
        background: #f8f9fa;
        border-radius: 8px;
        border: 1px solid #dee2e6;
        min-height: 50px;
    `;

    const startButton = document.createElement('button');
    startButton.textContent = 'Start Animation';
    startButton.style.cssText = `
        padding: 10px 20px;
        background: #2ecc71;
        color: white;
        border: none;
        border-radius: 4px;
        cursor: pointer;
        font-size: 16px;
    `;

    controlsContainer.appendChild(reconstructionDisplay);
    controlsContainer.appendChild(startButton);
    document.getElementById("graph").parentNode.insertBefore(controlsContainer, document.getElementById("graph"));

    const width = DEFAULT_WIDTH;
    const height = DEFAULT_HEIGHT;
    const margin = 40;

    // Clear existing graph
    d3.select("#graph").html("");

    const svg = d3.select("#graph")
        .append("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("viewBox", [0, 0, width, height])
        .attr("style", "background-color: #f8f9fa; border-radius: 8px;");

    // Add gradient definition for nodes
    const defs = svg.append("defs");
    const gradient = defs.append("radialGradient")
        .attr("id", "nodeGradient");
    
    gradient.append("stop")
        .attr("offset", "0%")
        .attr("stop-color", "#fff");
    
    gradient.append("stop")
        .attr("offset", "100%")
        .attr("stop-color", "#f0f0f0");

    // Arrow marker definition
    svg.append("defs").selectAll("marker")
        .data(["arrow"])
        .enter()
        .append("marker")
        .attr("id", d => d)
        .attr("viewBox", "0 -5 10 10")
        .attr("refX", 20)
        .attr("refY", 0)
        .attr("markerWidth", 8)
        .attr("markerHeight", 8)
        .attr("orient", "auto")
        .append("path")
        .attr("d", "M0,-5L10,0L0,5")
        .style("fill", "#666");

    const simulation = d3.forceSimulation(graphNodes)
        .force("link", d3.forceLink(graphLinks)
            .id(d => d.name)
            .distance(100))
        .force("charge", d3.forceManyBody().strength(-300))
        .force("center", d3.forceCenter(width / 2, height / 2))
        .force("collision", d3.forceCollide().radius(30))
        .force("x", d3.forceX(width / 2).strength(0.1))
        .force("y", d3.forceY(height / 2).strength(0.1));

    // Create paths for edges
    const path = svg.selectAll(".link")
        .data(graphLinks)
        .enter()
        .append("path")
        .attr("id", (d, i) => `link${i}`)
        .attr("class", "link")
        .attr("marker-end", "url(#arrow)")
        .style("stroke", "#666")
        .style("stroke-width", 1.5);

    // Add edge labels
    const edgeLabels = svg.selectAll(".edge-label")
        .data(graphLinks)
        .enter()
        .append("text")
        .attr("class", "edge-label")
        .append("textPath")
        .attr("xlink:href", (d, i) => `#link${i}`)
        .attr("startOffset", "50%")
        .style("text-anchor", "middle")
        .style("font-size", "12px")
        .style("fill", "#444")
        .text(d => d.transition || "");

    // Create node groups
    const node = svg.selectAll(".node")
        .data(graphNodes)
        .enter()
        .append("g")
        .attr("class", "node")
        .call(d3.drag()
            .on("start", dragstarted)
            .on("drag", dragged)
            .on("end", dragended));

    // Add circles to nodes
    node.append("circle")
        .attr("r", 10)
        .attr("class", "node-circle")
        .style("fill", d => d.name.includes('$') ? "#95a5a6" : "#e74c3c")
        .style("stroke", "#fff")
        .style("stroke-width", "2px");

    // Add node labels
    // First append a background rectangle for the text
    node.append("rect")
        .attr("class", "label-background")
        .attr("x", 12)
        .attr("y", -10)
        .attr("rx", 3)  // Rounded corners
        .attr("ry", 3)
        .attr("fill", "white")
        .attr("fill-opacity", 0.8)
        .attr("stroke", "#e5e5e5")
        .attr("stroke-width", 1);

    // Then append the text
    const labels = node.append("text")
        .attr("x", 12)
        .attr("y", ".31em")
        .style("font-family", "Arial, sans-serif")
        .style("font-size", "12px")
        .style("fill", "#333")
        .text(d => d.name);

    // Size the background rectangles to match the text
    node.selectAll(".label-background")
        .each(function(d) {
            const bbox = node.select("text").node().getBBox();
            d3.select(this)
                .attr("width", bbox.width + 6)   // Add padding
                .attr("height", bbox.height + 4)  // Add padding
                .attr("y", bbox.y - 2);          // Center vertically
        });

    function dragstarted(d) {
        if (!d3.event.active) simulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
    }

    function dragged(d) {
        d.fx = clamp(d3.event.x, margin, width - margin);
        d.fy = clamp(d3.event.y, margin, height - margin);
    }

    function dragended(d) {
        if (!d3.event.active) simulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
    }

    simulation.on("tick", () => {
        graphNodes.forEach(d => {
            d.x = clamp(d.x, margin, width - margin);
            d.y = clamp(d.y, margin, height - margin);
        });

        node.attr("transform", d => `translate(${d.x},${d.y})`);
        
        path.attr("d", d => {
            const dx = d.target.x - d.source.x;
            const dy = d.target.y - d.source.y;
            return `M${d.source.x},${d.source.y}L${d.target.x},${d.target.y}`;
        });
    });

    let animationTimeout;
    function startAnimation() {
        // Clear any existing animation
        clearTimeout(animationTimeout);
        
        // Reset display
        reconstructionDisplay.textContent = "Original string: " + originalString;
        
        // Reset visual states
        d3.selectAll(".node-circle").style("fill", d => d.name.includes('$') ? "#95a5a6" : "#e74c3c");
        d3.selectAll(".link").style("stroke", "#666").style("stroke-width", "1.5px");
        
        // Find path that matches original string
        const eulerianPath = findPathToMatch(graphNodes, graphLinks, originalString, k);
        animatePath(eulerianPath, node, path, reconstructionDisplay, k, originalString);
    }

    // Attach click handler to start button
    startButton.onclick = startAnimation;
}

function animatePath(path, nodes, edges, display, k, originalString) {
    let reconstructedString = originalString.substring(0, k - 1); // Start with the first k-1 characters
    let currentIndex = 0;

    display.innerHTML = `Original string: ${originalString}<br>Reconstructed: ${reconstructedString}`;

    function animate() {
        if (currentIndex >= path.length) {
            // Final validation and padding removal
            if (reconstructedString.endsWith('$')) {
                reconstructedString = reconstructedString.slice(0, -1);
            }
            display.innerHTML = `Original string: ${originalString}<br>Reconstructed:&nbsp&nbsp&nbsp${reconstructedString}`;
            return;
        }

        const step = path[currentIndex];

        // Highlight current node
        nodes.filter(n => n.name === step.from)
            .select("circle")
            .style("fill", "#2ecc71")
            .style("filter", "brightness(1.2)");

        // Highlight edge
        edges.filter(e => e.source.name === step.from && e.target.name === step.to)
            .style("stroke", "#2ecc71")
            .style("stroke-width", "3px");

        // Append the transition character
        reconstructedString += step.transition;
        display.innerHTML = `Original string: ${originalString}<br>Reconstructed:&nbsp&nbsp&nbsp${reconstructedString}`;

        currentIndex++;
        setTimeout(animate, 500);
    }

    animate();
}


function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
}