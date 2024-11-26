// Default dimensions
const DEFAULT_WIDTH = 800;
const DEFAULT_HEIGHT = 600;

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

function buildDeBruijnGraph(sequence, k) {
    const kmers = [];
    const links = [];
    const nodes = {};
    
    // Add proper padding
    const startPad = '$'.repeat(k-1);
    const paddedSequence = startPad + sequence + '$';
    
    // Process kmers
    for (let i = 0; i <= paddedSequence.length - k; i++) {
        const kmer = paddedSequence.substring(i, i + k);
        const prefix = kmer.substring(0, k-1);
        const suffix = kmer.substring(1, k);
        const transition = kmer[k-1];
        
        // Add to kmers list
        kmers.push({
            number: kmers.length + 1,
            length: prefix === startPad ? 0 : 1, // Length is 0 for pure $ nodes
            vertex: prefix,
            word: transition
        });

        // Create link
        links.push({
            source: prefix,
            target: suffix,
            transition: transition
        });
    }

    // Convert links to node references
    links.forEach(function(link) {
        link.source = nodes[link.source] || (nodes[link.source] = {name: link.source});
        link.target = nodes[link.target] || (nodes[link.target] = {name: link.target});
    });

    return { nodes: d3.values(nodes), links, kmers };
}

function renderKmerTable(kmers) {
    let tableContainer = document.getElementById('kmer-table');
    if (!tableContainer) {
        tableContainer = document.createElement('div');
        tableContainer.id = 'kmer-table';
        tableContainer.style.position = 'absolute';
        tableContainer.style.right = '20px';
        tableContainer.style.top = '20px';
        tableContainer.style.backgroundColor = '#f5f5f5';
        tableContainer.style.padding = '10px';
        tableContainer.style.border = '1px solid #ccc';
        tableContainer.style.borderRadius = '4px';
        document.querySelector('.container').appendChild(tableContainer);
    }

    const table = document.createElement('table');
    table.style.borderCollapse = 'collapse';
    table.style.width = '200px';
    table.style.fontSize = '14px';
    
    // Create header
    const header = table.createTHead();
    const headerRow = header.insertRow();
    ['#', 'L', 'Ver', 'W'].forEach(text => {
        const th = document.createElement('th');
        th.textContent = text;
        th.style.border = '1px solid #ccc';
        th.style.padding = '5px 10px';
        th.style.backgroundColor = '#e0e0e0';
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
            td.style.border = '1px solid #ccc';
            td.style.padding = '5px 10px';
            td.style.textAlign = 'center';
            row.appendChild(td);
        });
    });

    tableContainer.innerHTML = '';
    tableContainer.appendChild(table);
}

function clamp(x, min, max) {
    return Math.max(min, Math.min(max, x));
}

function renderGraph(graphNodes, graphLinks) {
    const width = document.getElementById('viewer')?.clientWidth || DEFAULT_WIDTH;
    const height = document.getElementById('viewer')?.clientHeight || DEFAULT_HEIGHT;
    const margin = 10;

    // Clear existing graph
    d3.select("#graph").html("");

    // Arrow marker definition
    const svg = d3.select("#graph")
        .append("svg")
        .attr("width", width)
        .attr("height", height);

    svg.append("defs").selectAll("marker")
        .data(["arrow"])
        .enter()
        .append("marker")
        .attr("id", d => d)
        .attr("viewBox", "0 -5 10 10")
        .attr("refX", 15) // Adjust based on circle radius
        .attr("refY", 0)
        .attr("markerWidth", 6)
        .attr("markerHeight", 6)
        .attr("orient", "auto")
        .append("path")
        .attr("d", "M0,-5L10,0L0,5") // Triangle path
        .style("fill", "#000");

    // Create force simulation
    const simulation = d3.forceSimulation(graphNodes)
        .force("link", d3.forceLink(graphLinks)
            .id(d => d.name)
            .distance(60))
        .force("charge", d3.forceManyBody().strength(-200))
        .force("center", d3.forceCenter(width / 2, height / 2));

    // Create paths for edges
    const path = svg.selectAll(".link")
        .data(graphLinks)
        .enter()
        .append("path")
        .attr("id", (d, i) => `link${i}`) // Unique ID for each link
        .attr("class", d => (d.target.name.includes('$') ? "link auxiliary" : "link"))
        .attr("marker-end", "url(#arrow)") // Add arrow to end of the path
        .style("stroke", "#000")
        .style("stroke-width", 1)
        .style("fill", "none");

    // Add edge labels
    svg.selectAll(".edge-label")
        .data(graphLinks)
        .enter()
        .append("text")
        .attr("class", "edge-label")
        .append("textPath")
        .attr("xlink:href", (d, i) => `#link${i}`) // Link to the path ID
        .attr("startOffset", "50%") // Position at the middle of the path
        .style("text-anchor", "middle")
        .style("font-size", "10px")
        .style("fill", "#000")
        .text(d => d.transition || ""); // Text for edge (transition character)

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
        .attr("r", 6)
        .style("fill", d => {
            if (d.name === "$".repeat(d.name.length)) {
                return "black";
            } else if (d.name.includes('$')) {
                return d.name[0] !== '$' ? "white" : "gray";
            } else {
                return "red";
            }
        })
        .style("stroke", d => (d.name.includes('$') && d.name[0] !== '$' ? "black" : "none"))
        .style("stroke-width", "0.5px");

    // Add node labels
    node.append("text")
        .attr("x", 8)
        .attr("y", ".31em")
        .text(d => d.name);

    // Create a drag handler and append it to the node object instead
    var drag_handler = d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended);
    drag_handler(node);

    function dragstarted(d) {
        if (!d3.event.active) simulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
    }

    function dragged(d) {
        d.fx = d3.event.x;
        d.fy = d3.event.y;
    }

    function dragended(d) {
        if (!d3.event.active) simulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
    }

    simulation.on("tick", () => {
        // Update node positions
        node.attr("transform", d => `translate(${d.x},${d.y})`);

        // Update edge positions
        path.attr("d", d => {
            const dx = d.target.x - d.source.x;
            const dy = d.target.y - d.source.y;
            const dr = Math.sqrt(dx * dx + dy * dy);
            return `M${d.source.x},${d.source.y}L${d.target.x},${d.target.y}`;
        });
    });

    function clamp(value, min, max) {
        return Math.max(min, Math.min(max, value));
    }

    // Add CSS styles
    const style = document.createElement('style');
    style.textContent = `
        .link {
            stroke: #000;
            stroke-width: 1px;
        }
        .link.auxiliary {
            stroke-dasharray: 2,2;
        }
        .node circle {
            fill-opacity: 1;
        }
        .node text {
            font: 10px sans-serif;
        }
        .edge-label {
            font: 10px sans-serif;
        }
    `;
    document.head.appendChild(style);
}



// // Add window resize handler
// window.addEventListener('resize', function() {
//     const viewer = document.getElementById('viewer');
//     if (viewer && viewer.querySelector('svg')) {
//         const width = viewer.clientWidth;
//         const height = viewer.clientHeight;
//         d3.select("#graph svg")
//             .attr("width", width)
//             .attr("height", height)
//             .attr("viewBox", [0, 0, width, height]);
//     }
// });