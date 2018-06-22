function computeAccelerations() {

    //
    var dx, dy, dx2, dy2, r, rSquared, rSquaredInv, attract, repel, fOverR, fx, fy;

    // forceCutoff is the distance that you stop trying to compute force
    var forceCutoff = 3.0;						// distance beyond which we set force=0

    // Square of forceCutoff
    var forceCutoff2 = forceCutoff * forceCutoff;

    // PE at cutoff
    var pEatCutoff = 4 * (Math.pow(forceCutoff, -12) - Math.pow(forceCutoff, -6));

    // Gravity
    var g = Number(getGravity());

    // Wall is a spring
    var wallStiffness = 50;						// spring constant for bouncing off walls

    // Variable to hold the force imparted by wall
    var wallForce = 0.0;

   // Starting PE 
    potentialE = 0.0;

    // first check for bounces off walls:
    for (var i = 0; i < N; i++) {
        if (x[i] < 0.5) {
            ax[i] = wallStiffness * (0.5 - x[i]);
            wallForce += ax[i];
            potentialE += 0.5 * wallStiffness * (0.5 - x[i]) * (0.5 - x[i]);
        } else
            if (x[i] > (boxWidth - 0.5)) {
                ax[i] = wallStiffness * (boxWidth - 0.5 - x[i]);
                wallForce -= ax[i];
                potentialE += 0.5 * wallStiffness * (boxWidth - 0.5 - x[i]) * (boxWidth - 0.5 - x[i]);
            } else
                ax[i] = 0.0;
        if (y[i] < 0.5) {
            ay[i] = (wallStiffness * (0.5 - y[i]));
            wallForce += ay[i];
            potentialE += 0.5 * wallStiffness * (0.5 - y[i]) * (0.5 - y[i]);
        } else
            if (y[i] > (boxWidth - 0.5)) {
                ay[i] = (wallStiffness * (boxWidth - 0.5 - y[i]));
                wallForce -= ay[i];
                potentialE += 0.5 * wallStiffness * (boxWidth - 0.5 - y[i]) * (boxWidth - 0.5 - y[i]);
            } else
                ay[i] = 0;
        ay[i] -= g;				// add gravity if any
    }
    pressure = wallForce / (4 * boxWidth);	// instantaneous pressure

    // now compute interaction forces (Lennard-Jones potential):
    // (this is where we spend most of our computation time, so try to optimize)
    
    // Looks like they've set some computation threshholds so that if the force cutoff is too small, it doesn't bother
    if ((N < 100) || (boxWidth < 4 * forceCutoff) /* || !cellListCheck.checked */) {

        // Iterate over all particles with an i and j index
        for (var i = 0; i < N; i++) {			// simple double-loop over atoms for small system
            for (var j = 0; j < i; j++) {
                
                // X-distance between i and jth particle, r_x
                dx = x[i] - x[j];

                // Square of the x-distance, r_x^2
                dx2 = dx * dx;

                // If they are close enough to one another, namely within force cutoff
                if (dx2 < forceCutoff2) {  // make sure they're close enough to bother

                    // Y-distance between i and jth particle, r_y
                    dy = y[i] - y[j];

                    // Square of the y-distance, r_y^2
                    dy2 = dy * dy;

                    // Is the y-distance also close enough, within force cutoff
                    if (dy2 < forceCutoff2) {

                        // Calculating the distance r^2
                        rSquared = dx2 + dy2;

                        // Is r^2 within force cutoff?
                        if (rSquared < forceCutoff2) {

                            // Calculates 1/r^2
                            rSquaredInv = 1.0 / rSquared;

                            // Attraction term = (1/r^2)^3 = 1/r^6
                            attract = rSquaredInv * rSquaredInv * rSquaredInv;

                            // Repulsion term = 1/r^6 * 1/r^6
                            repel = attract * attract;

                            // delta V relative to the cutoff value 
                            potentialE += (4.0 * (repel - attract)) - pEatCutoff;

                            // F/r is the force per unit distance
                            fOverR = 24.0 * ((2.0 * repel) - attract) * rSquaredInv;
                            // 24.0 * ((2.0 * 1/r^12) - 1/r^6) * 1/r^2 = 24.0 * ((2.0 * 1/r^13) - 1/r^7) * 1/r

                            fx = fOverR * dx;
                            fy = fOverR * dy;
                            ax[i] += fx;  // add this force on to i's acceleration (m = 1)
                            ay[i] += fy;
                            ax[j] -= fx;  // Newton's 3rd law
                            ay[j] -= fy;
                        }
                    }
                }
            }
        }
    } else {	// tricky O(N) cell-based approach for large system
        var nCells, cellWidth, xCell, yCell, thisCell, neighborCell, xNeighborCell, yNeighborCell;
        var neighborOffset = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 1, y: 1 }, { x: 0, y: 1 }, { x: -1, y: 1 }];	// here, E, NE, N, and NW
        nCells = Math.floor(boxWidth / forceCutoff);		// number of cells in a row
        cellWidth = boxWidth / nCells;
        var listHeader = new Array(nCells * nCells);			// linked list headers
        for (var cell = 0; cell < nCells * nCells; cell++) listHeader[cell] = -1;		// set all cells to empty
        var linkedList = new Array(N);		// element i will point to next atom in same cell
        for (var i = 0; i < N; i++) {			// this loop assembles the linked list of atoms by cell
            xCell = Math.floor(x[i] / cellWidth);		// figure out which cell the atom is in
            if (xCell < 0) xCell = 0;
            if (xCell >= nCells) xCell = nCells - 1;
            yCell = Math.floor(y[i] / cellWidth);
            if (yCell < 0) yCell = 0;
            if (yCell >= nCells) yCell = nCells - 1;
            var cellHeaderIndex = xCell + nCells * yCell;		// flatten 2D structure into 1D array
            linkedList[i] = listHeader[cellHeaderIndex];	// this atom now points where the header used to
            listHeader[cellHeaderIndex] = i;				// header now points to his atom
        }	// linked list is now complete
        for (xCell = 0; xCell < nCells; xCell++) {				// loop over cells
            for (yCell = 0; yCell < nCells; yCell++) {
                thisCell = xCell + nCells * yCell;		// index of this cell in header list
                for (var neighborIndex = 0; neighborIndex < 5; neighborIndex++) {	// loop over neighboring cells
                    xNeighborCell = xCell + neighborOffset[neighborIndex].x;
                    if ((xNeighborCell < 0) || (xNeighborCell >= nCells)) continue;	// some neighbors don't actually exist
                    yNeighborCell = yCell + neighborOffset[neighborIndex].y;
                    if (yNeighborCell >= nCells) continue;
                    neighborCell = xNeighborCell + nCells * yNeighborCell;	// index of neighbor cell in header list
                    var i = listHeader[thisCell];
                    while (i > -1) {	// loop over atoms in this cell
                        var j = listHeader[neighborCell];
                        if (neighborCell == thisCell) j = linkedList[i];	// be sure not to count atoms in this cell twice
                        while (j > -1) {	// loop over atoms in neighbor cell
                            dx = x[i] - x[j];
                            dx2 = dx * dx;
                            if (dx2 < forceCutoff2) {  // make sure they're close enough to bother
                                dy = y[i] - y[j];
                                dy2 = dy * dy;
                                if (dy2 < forceCutoff2) {
                                    rSquared = dx2 + dy2;
                                    if (rSquared < forceCutoff2) {
                                        rSquaredInv = 1.0 / rSquared;
                                        attract = rSquaredInv * rSquaredInv * rSquaredInv;
                                        repel = attract * attract;
                                        potentialE += (4.0 * (repel - attract)) - pEatCutoff;
                                        fOverR = 24.0 * ((2.0 * repel) - attract) * rSquaredInv;
                                        fx = fOverR * dx;
                                        fy = fOverR * dy;
                                        ax[i] += fx;  // add this force on to i's acceleration (m = 1)
                                        ay[i] += fy;
                                        ax[j] -= fx;  // Newton's 3rd law
                                        ay[j] -= fy;
                                    }
                                }
                            }
                            j = linkedList[j];
                        } // end of loop over j
                        i = linkedList[i];
                    } // end of loop over i
                } // end of loop over neighborIndex
            } // end of loop over yCell
        } // end of loop over xCell
    } // end if (and end of L-J force computation)

    // add elastic forces between bonded atoms:
    var bondStrength = 100;	// spring constant (vastly less than actual covalent bonds!)
    for (var i = 0; i < bondCount * 2; i += 2) {
        var i1 = bondList[i];
        var i2 = bondList[i + 1];
        dx = x[i1] - x[i2];
        dy = y[i1] - y[i2];
        r = Math.sqrt(dx * dx + dy * dy);
        var rOffset = r - 1.122462;		// offset by L-J equilibrium position
        potentialE += 0.5 * bondStrength * rOffset * rOffset;
        fx = bondStrength * rOffset * dx / r;
        fy = bondStrength * rOffset * dy / r;
        ax[i1] -= fx;
        ay[i1] -= fy;
        ax[i2] += fx;
        ay[i2] += fy;
    }

    // fixed atoms don't accelerate:
    for (var i = 0; i < fixedCount; i++) {
        ax[fixedList[i]] = 0;
        ay[fixedList[i]] = 0;
    }

    // if we're pulling on an atom it feels an elastic force toward mouse location:
    if (dragging) {
        var pullStrength = 1.0;			// spring constant
        dx = mouseX - x[clickedAtom];
        dy = mouseY - y[clickedAtom];
        ax[clickedAtom] += pullStrength * dx;
        ay[clickedAtom] += pullStrength * dy;
    }
}