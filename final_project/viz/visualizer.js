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
    renderGraph(nodes, links);
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
            transition: transition
        });
    }

    links.forEach(function(link) {
        link.source = nodes[link.source] || (nodes[link.source] = {name: link.source});
        link.target = nodes[link.target] || (nodes[link.target] = {name: link.target});
    });

    return { nodes: d3.values(nodes), links, kmers };
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

function renderGraph(graphNodes, graphLinks) {
    const width = DEFAULT_WIDTH;
    const height = DEFAULT_HEIGHT;
    const margin = 40; // Increased margin for better spacing

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

    // Arrow marker definition with improved styling
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

    // Create force simulation with bounds
    const simulation = d3.forceSimulation(graphNodes)
        .force("link", d3.forceLink(graphLinks)
            .id(d => d.name)
            .distance(100))
        .force("charge", d3.forceManyBody().strength(-300))
        .force("center", d3.forceCenter(width / 2, height / 2))
        .force("collision", d3.forceCollide().radius(30))
        .force("x", d3.forceX(width / 2).strength(0.1))
        .force("y", d3.forceY(height / 2).strength(0.1));

    // Create paths for edges with improved styling
    const path = svg.selectAll(".link")
        .data(graphLinks)
        .enter()
        .append("path")
        .attr("id", (d, i) => `link${i}`)
        .attr("class", d => (d.target.name.includes('$') ? "link auxiliary" : "link"))
        .attr("marker-end", "url(#arrow)")
        .style("stroke", "#666")
        .style("stroke-width", 1.5);

    // Add edge labels with improved styling
    svg.selectAll(".edge-label")
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

    // Create node groups with improved styling
    const node = svg.selectAll(".node")
        .data(graphNodes)
        .enter()
        .append("g")
        .attr("class", "node")
        .call(d3.drag()
            .on("start", dragstarted)
            .on("drag", dragged)
            .on("end", dragended));

    // Add circles to nodes with improved styling
    node.append("circle")
        .attr("r", 10)
        .style("fill", d => {
            if (d.name === "$".repeat(d.name.length)) {
                return "#2c3e50";
            } else if (d.name.includes('$')) {
                return d.name[0] !== '$' ? "url(#nodeGradient)" : "#95a5a6";
            } else {
                return "#e74c3c";
            }
        })
        .style("stroke", "#fff")
        .style("stroke-width", "2px")
        .style("filter", "drop-shadow(0 2px 2px rgba(0,0,0,0.1))");

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
        // Keep nodes within bounds
        graphNodes.forEach(d => {
            d.x = clamp(d.x, margin, width - margin);
            d.y = clamp(d.y, margin, height - margin);
        });

        // Update node positions
        node.attr("transform", d => `translate(${d.x},${d.y})`);

        // Update edge positions with curved paths
        path.attr("d", d => {
            const dx = d.target.x - d.source.x;
            const dy = d.target.y - d.source.y;
            return `M${d.source.x},${d.source.y}L${d.target.x},${d.target.y}`;
        });
    });

    // Add CSS styles
    const style = document.createElement('style');
    style.textContent = `
        .link {
            stroke: #666;
            stroke-width: 1.5px;
        }
        .link.auxiliary {
            stroke-dasharray: 4,4;
            stroke-opacity: 0.6;
        }
        .node circle {
            fill-opacity: 1;
            transition: all 0.3s ease;
        }
        .node:hover circle {
            filter: brightness(1.1);
        }
        .node text {
            font-family: Arial, sans-serif;
            font-size: 12px;
            pointer-events: none;
        }
        .edge-label {
            font-family: Arial, sans-serif;
            pointer-events: none;
        }
    `;
    document.head.appendChild(style);
}

function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
}

// Add window resize handler
window.addEventListener('resize', function() {
    const viewer = document.getElementById('viewer');
    if (viewer && viewer.querySelector('svg')) {
        const width = Math.max(DEFAULT_WIDTH, viewer.clientWidth);
        const height = Math.max(DEFAULT_HEIGHT, viewer.clientHeight);
        d3.select("#graph svg")
            .attr("width", width)
            .attr("height", height)
            .attr("viewBox", [0, 0, width, height]);
    }
});