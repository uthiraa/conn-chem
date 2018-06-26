var moleculeArray = [];
var boundary;
var piston;

var systemKineticEnergy;
var uniqueForcePairs;

var imageDir = "img/svg/";

var numberOfParticles = 20;

function preload() {

}

function setup() {
    // Creates the drawing canvas
    createCanvas(800, 600);

    // Sets all x,y coordinates to the center of images
    imageMode(CENTER);

    // Adjust framerate to slow animation, useful for debugging
    // frameRate(5);

    // Initiates a new b2 world, scaling factor = 30, gravity vector
    b2newWorld(30, createVector(0, 9.8));

    // Creates the boundary for a closed reaction container
    boundary = new Boundary("CLOSED");

    // Creates the piston boundary for the reaction container
    piston = new Boundary("PISTON", createVector(width / 2, 0), createVector(width, 20))

    // Calculates how many rows of particles to draw, currently assumes 5 columns of particles
    var rowsOfParticles = Math.floor(numberOfParticles / 5)

    // Increments the number of rows if there isn't exactly a multiple of 5 particles
    if (numberOfParticles % 5 > 0) {
        rowsOfParticles++;
    }

    // Creates molecule collection in a grid with random velocity vectors
    for (var i = 0; i < rowsOfParticles; i++) {
        for (var j = 0; j < 5; j++) {
            moleculeArray[5 * i + j] = new Particle(161, createVector(width / 10 + (j * (width / 5)), height / 10 + i * (height / 10)), p5.Vector.random2D().mult(random(0, 1)), 5 * i + j);
        }
    }
}

function draw() {
    // The background must be redrawn with each draw loop to avoid molecules leaving permanent traces
    background(204, 204, 204);

    // Resets the KE of the system each draw loop
    systemKineticEnergy = 0;

    // Physics calculations must be updated each draw by calling b2Update
    b2Update();
    b2Draw(false);

    // Calculate intermolecular forces and calculate net force vector
    // This is currently not working properly
    for (var i = 0; i < numberOfParticles; i++) {
        for (var j = 0; j < numberOfParticles; j++) {
            if (moleculeArray[i].id == moleculeArray[j].id || moleculeArray[i].indexedBy(moleculeArray[j])) {
                continue;
            } else {
                moleculeArray[i].indexes(moleculeArray[j]);

                // Apply the force from the vdW interactions (to be implemented)
                var ljVector = moleculeArray[i].calculateLJForce(moleculeArray[j], 1e2);

                // Newton's 3rd Law
                moleculeArray[i].addForceToNetForce(ljVector);
                moleculeArray[j].addForceToNetForce(ljVector.mult(-1));

            }
        }
    }

    // Iterates over particles to update system KE & apply forces
    for (var i = 0; i < numberOfParticles; i++) {
        // Calculates the KE of each particle
        moleculeArray[i].updateKineticEnergy();

        // KE for each particle is added to the system KE for display
        systemKineticEnergy += moleculeArray[i].getKineticEnergy();

        // Applies intermolecular force calculations
        moleculeArray[i].applyNetForce();
        // console.log(moleculeArray[i].getNetForce());
        // moleculeArray[i].showNetForce();

        // Resets the net force to 0 and clears out force indices each draw loop
        moleculeArray[i].clearForceIndices();
        moleculeArray[i].clearNetForce();
    }

    // Displays a readout of the system KE
    text(systemKineticEnergy, width / 2, height - 10);
}

// Boundary class to create common boundary types
class Boundary {

    constructor(type, position, dimensions) {
        // Creates a closed container to bound the entire system
        if (type.toUpperCase() == "CLOSED") {
            this.bottom = new b2Body('box', false, createVector(width / 2, height), createVector(width, 1));
            this.left = new b2Body('box', false, createVector(0, height / 2), createVector(1, height));
            this.right = new b2Body('box', false, createVector(width, height / 2), createVector(1, height));
            this.top = new b2Body('box', false, createVector(width / 2, 0), createVector(width, 1));
        }

        // Creates a horizontal boundary for the piston 
        if (type.toUpperCase() == "PISTON") {
            this.piston = new b2Body('box', false, createVector(position.x, position.y), createVector(dimensions.x, dimensions.y));
            this.piston.display(this.drawBoundary, 0);
        }
    }

    drawBoundary(body, fixture, position) {
        fill(150, 150, 150);
        rect(position.x, position.y, body.wh(0).x, body.wh(0).y)
        // b2Display(body, fixture, position);
    }
}

// Particle class to handle each particle and its properties
class Particle {

    constructor(key, position, velocity, id) {
        // This is the value in the "ID" field of the database variable
        this.databaseKey = key;

        // This gives each molecule a unique index to track it within a set of molecules
        this.id = id;

        // Particle ID is database key, used to pull name
        this.name = database[this.databaseKey - 1].name;

        // Sets the position of the particle
        this.position = createVector(position.x, position.y);

        // Sets the velocity of the particle
        this.velocity = createVector(velocity.x, velocity.y)

        // Initializes the particle mass
        this.mass = 1;

        // When forces are applied, we need to indicate which particle force pairs have already been calculated
        // Avoids duplicate calculation of intermolecular forces
        this.forceIndices = [];

        // Net force vector on the particle
        this.netForce = createVector(0, 0);

        // Sets the KE of the particle
        this.kineticEnergy = 0;

        // Sets the image reference for using the ID number of the particle
        this.imageUrl = imageDir + database[this.databaseKey - 1].file;
        this.imageObject = loadImage(this.imageUrl, (result) => {

            // We have to invoke the callback function since loading images are done asynchronously
            // The result, the image object, is then passed to the drawParticle method
            this.drawParticle(result);

        });

    }

    /**
     * Calculate properties of particle
     */

    updateKineticEnergy() {

        // Check first that this.body is defined
        if (!(typeof this.body == 'undefined')) {

            // Then the kinetic energy is calculated
            this.kineticEnergy = this.getTranslationalKineticEnergy() + this.getRotationalKineticEnergy();
        }
    }

    // Returns distance vector starting at this particle and pointing toward the particle in the argument
    vectorToParticle(particle) {
        return p5.Vector.sub(particle.getPosition(), this.getPosition());
    }

    // Returns the magnitude of the distance vector between two particles
    distanceToParticle(particle) {
        return this.vectorToParticle(particle).mag();
    }

    /**
     * Calculate forces on the particle
     */

    // Adds a vector of a force affecting a particle 
    addForceToNetForce(force) {
        this.netForce.add(force);
    }

    // Applies a force onto a particle using the b2 library
    applyNetForce() {
        if (!(typeof this.body == 'undefined')) {
            var magnitude = this.netForce.mag();
            this.body.applyForce(createVector(this.netForce.x / magnitude, this.netForce.y / magnitude), magnitude);
        }
    }

    // Returns whether the argument particle applied a force to this particle
    indexedBy(particle) {
        return this.forceIndices.includes(particle.id);
    };

    // Adds the indices of this and the particle argument to the force index array for book keeping
    indexes(particle) {
        this.forceIndices.push(particle.id);
        particle.forceIndices.push(this.id);
    }

    // Clears out the book keeping force pairs after each loop
    clearForceIndices() {
        this.forceIndices = [];
    }

    // Sets the net force back to zero after each iteration
    clearNetForce() {
        this.netForce.x = 0;
        this.netForce.y = 0;
    }

    // Calculates the van der Waals force on a particle using a Lennard-Jones potential
    // The Lennard-Jones potential is a mathematical simplification of the potential energy caused by van der Waals interactions
    calculateLJForce(particle, epsilon, forceCutoff) {
        var epsilonValue, forceCutoffValue, forceToReturn;
        var distanceBetweenParticles = this.distanceToParticle(particle);
        var vectorBetweenParticles = this.vectorToParticle(particle);

        // Sets the depth of the potential well
        if (typeof epsilon == "undefined") {
            epsilonValue = 1.0;
        } else {
            epsilonValue = epsilon;
        }

        // Sets a maximum force, attempt to prevent force blow-up at small delta r
        if (typeof forceCutoff == "undefined") {
            forceCutoffValue = 1.0;
        } else {
            forceCutoffValue = forceCutoff;
        }

        // Determines whether the force should be applied
        if (distanceBetweenParticles > 10.0 * this.getRadius()) {
            forceToReturn = 0;
        } else {
            // Technique #1: Using absolute particle distance
            // 24.0 * ((2.0 * 1/r^14) - 1/r^8) = 24.0 * ((2.0 * 1/r^13) - 1/r^7) * 1/r = F(r) * 1/r
            // var fOverR = 24.0 * epsilonValue * ((2.0 * Math.pow(distanceBetweenParticles, -14)) - Math.pow(distanceBetweenParticles, -8));

            // Technique #2: Using ratio of sprite radius (max of length, width) to particle distance
            var fOverR = 12 * epsilonValue * (Math.pow((this.getRadius() / distanceBetweenParticles),13) - Math.pow((this.getRadius() / distanceBetweenParticles),7)) * Math.pow(distanceBetweenParticles, -1);

            forceToReturn = Math.min(fOverR * distanceBetweenParticles, forceCutoffValue);            
        }
        // Negative forces are attractive by definition, so we multiply by -1 to ensure the force vector has the same direction as r
        return vectorBetweenParticles.normalize().mult(-forceToReturn);
    }

    /**
     * Get functions
     */

    // Returns a P5 Vector containing the particle velocity
    getVelocity() {
        return (!(typeof this.body == 'undefined')) ? this.body.velocity : this.velocity;
    }

    // Returns the angular velocity of a particle
    getAngularVelocity() {
        return (!(typeof this.body == 'undefined')) ? this.body.angularVelocity : 0;
    }

    // Returns the position in pixels of the particle
    getPosition() {
        return (!(typeof this.body == 'undefined')) ? this.body.xy : this.position;
    }

    // Returns the mass of the body
    getMass() {
        return (!(typeof this.body == 'undefined')) ? this.body.mass : this.mass;
    }

    // Returns the kinetic energy (rotational + translational) of the particle
    getKineticEnergy() {
        return this.kineticEnergy;
    }

    // Returns the translational kinetic energy of the particle
    getTranslationalKineticEnergy() {
        return (0.5 * this.body.mass * Math.pow(this.getVelocity().mag(), 2));
    }

    // Returns the rotational kinetic energy of the particle
    getRotationalKineticEnergy() {
        return ((1 / 24) * this.getMass() * Math.pow(this.getAngularVelocity(), 2) * (Math.pow(this.imageObject.width, 2) + Math.pow(this.imageObject.height, 2)));
    }

    // Takes the maximum dimension of the molecule sprite as an approximate radial extent of the molecule
    getRadius() {
        return (!(typeof this.body == 'undefined')) ? Math.max(this.body.wh(0).x / 2.0, this.body.wh(0).y / 2.0) : Math.max(this.imageObject.width / 2.0, this.imageObject.height / 2.0);
    }

    // Returns the net force vector (intermolecular forces)
    getNetForce() {
        return this.netForce;
    }

    /**
     * Set functions
     */

    // Updates the vx,vy based on the force field
    setVelocity(velocity) {
        this.velocity.x = velocity.x;
        this.velocity.y = velocity.y;
    }

    // Sets the x,y position based on the arguments
    setPosition(position) {
        this.position.x = position.x;
        this.position.y = position.y;
    }

    /*
    * Rendering functions
    */

    // Renders the image glyph to canvas using P5.js Image object
    drawParticle(imageObject) {
        // Creates b2 Body to interact with b2 world
        // This is called asychronously after the image is loaded
        this.body = new b2Body('type', true, createVector(this.position.x, this.position.y), createVector(this.imageObject.width, this.imageObject.height), 1.0, 0, 1.0);
        this.body.image(imageObject, 0);

        // Set particle velocity
        this.body.applyForce(createVector(this.velocity.x, this.velocity.y), 100);

        // Sets position on the class variable
        this.setPosition(this.body.xy);
    }

    // Displays a velocity vector representation for each particle
    showVelocity() {
        line(this.position.x, this.position.y, (this.position.x + 3 * this.velocity.x), (this.position.y + 3 * this.velocity.y));
    }

    // Displays a hovering bookeeping index for each particle
    showLabel() {
        text(this.id, this.position.x + 25, this.position.y + 25);
    }

    // Shows a line representation of the net force vector
    showNetForce(visualScale) {
        var scale = (typeof visualScale == "undefined") ? 3 : visualScale;

        line(this.getPosition().x, this.getPosition().y, (this.getPosition().x + scale * this.netForce.x), (this.getPosition().y + scale * this.netForce.y));
    }

    /**
     * Logging functions
    */

    // Writes particle ID, position, and velocity to the console
    log() {
        // console.log(this.id + "\t" + this.getPosition().x + "\t" + this.getPosition().y + "\t" + this.getVelocity().x + "\t" + this.getVelocity().y);
        if (typeof (this.body) == "undefined") {
            // Do nothing - this is necessary to wait until the object is created
        } else {
            console.log(this.body.velocity);
        }
    }

    // Provides a string representation of the Particle object
    toString() {
        return "This particle has the chemical identity '" + this.name + "'";
    }

};

// Array of Javascript objects containing names and file locations for different species
var database = [
    {
        "ID": 1,
        "file": "384100F6.png",
        "name": "384100F6",
        "nameSlug": "384100f6",
        "fileExtension": "png"
    },
    {
        "ID": 2,
        "file": "6161D219.png",
        "name": "6161D219",
        "nameSlug": "6161d219",
        "fileExtension": "png"
    },
    {
        "ID": 3,
        "file": "6161D21A.png",
        "name": "6161D21A",
        "nameSlug": "6161d21a",
        "fileExtension": "png"
    },
    {
        "ID": 4,
        "file": "Acetate.svg",
        "name": "Acetate",
        "nameSlug": "acetate",
        "fileExtension": "svg"
    },
    {
        "ID": 5,
        "file": "Acetic-Acid.svg",
        "name": "Acetic Acid",
        "nameSlug": "acetic-acid",
        "fileExtension": "svg"
    },
    {
        "ID": 6,
        "file": "Aluminum-Chloride.svg",
        "name": "Aluminum Chloride",
        "nameSlug": "aluminum-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 7,
        "file": "Aluminum-Ion.svg",
        "name": "Aluminum Ion",
        "nameSlug": "aluminum-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 8,
        "file": "Aluminum-Oxide.svg",
        "name": "Aluminum Oxide",
        "nameSlug": "aluminum-oxide",
        "fileExtension": "svg"
    },
    {
        "ID": 9,
        "file": "Aluminum.svg",
        "name": "Aluminum",
        "nameSlug": "aluminum",
        "fileExtension": "svg"
    },
    {
        "ID": 10,
        "file": "Ammonia-Dots.svg",
        "name": "Ammonia Dots",
        "nameSlug": "ammonia-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 11,
        "file": "Ammonia.svg",
        "name": "Ammonia",
        "nameSlug": "ammonia",
        "fileExtension": "svg"
    },
    {
        "ID": 12,
        "file": "Ammonium-Acetate.svg",
        "name": "Ammonium Acetate",
        "nameSlug": "ammonium-acetate",
        "fileExtension": "svg"
    },
    {
        "ID": 13,
        "file": "Ammonium-Chloride.svg",
        "name": "Ammonium Chloride",
        "nameSlug": "ammonium-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 14,
        "file": "Ammonium-Dots.svg",
        "name": "Ammonium Dots",
        "nameSlug": "ammonium-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 15,
        "file": "Ammonium-Nitrate.svg",
        "name": "Ammonium Nitrate",
        "nameSlug": "ammonium-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 16,
        "file": "Ammonium.svg",
        "name": "Ammonium",
        "nameSlug": "ammonium",
        "fileExtension": "svg"
    },
    {
        "ID": 17,
        "file": "Barium-144.svg",
        "name": "Barium 144",
        "nameSlug": "barium-144",
        "fileExtension": "svg"
    },
    {
        "ID": 18,
        "file": "Bicarbonate.svg",
        "name": "Bicarbonate",
        "nameSlug": "bicarbonate",
        "fileExtension": "svg"
    },
    {
        "ID": 19,
        "file": "Blue.svg",
        "name": "Blue",
        "nameSlug": "blue",
        "fileExtension": "svg"
    },
    {
        "ID": 20,
        "file": "Boron-Ion.svg",
        "name": "Boron Ion",
        "nameSlug": "boron-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 21,
        "file": "Boron.svg",
        "name": "Boron",
        "nameSlug": "boron",
        "fileExtension": "svg"
    },
    {
        "ID": 22,
        "file": "Boron-Tetrachloride-Dots.svg",
        "name": "Boron Tetrachloride Dots",
        "nameSlug": "boron-tetrachloride-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 23,
        "file": "Boron-Tetrachloride.svg",
        "name": "Boron Tetrachloride",
        "nameSlug": "boron-tetrachloride",
        "fileExtension": "svg"
    },
    {
        "ID": 24,
        "file": "Boron-Trichloride-Dots.svg",
        "name": "Boron Trichloride Dots",
        "nameSlug": "boron-trichloride-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 25,
        "file": "Boron-Trichloride.svg",
        "name": "Boron Trichloride",
        "nameSlug": "boron-trichloride",
        "fileExtension": "svg"
    },
    {
        "ID": 26,
        "file": "Bromide-Dots.svg",
        "name": "Bromide Dots",
        "nameSlug": "bromide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 27,
        "file": "Bromide.svg",
        "name": "Bromide",
        "nameSlug": "bromide",
        "fileExtension": "svg"
    },
    {
        "ID": 28,
        "file": "Bromine-Atom.svg",
        "name": "Bromine Atom",
        "nameSlug": "bromine-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 29,
        "file": "Bromine-Ion.svg",
        "name": "Bromine Ion",
        "nameSlug": "bromine-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 30,
        "file": "Bromine.svg",
        "name": "Bromine",
        "nameSlug": "bromine",
        "fileExtension": "svg"
    },
    {
        "ID": 31,
        "file": "Butane.svg",
        "name": "Butane",
        "nameSlug": "butane",
        "fileExtension": "svg"
    },
    {
        "ID": 32,
        "file": "Butene.svg",
        "name": "Butene",
        "nameSlug": "butene",
        "fileExtension": "svg"
    },
    {
        "ID": 33,
        "file": "Calcium-Chloride.svg",
        "name": "Calcium Chloride",
        "nameSlug": "calcium-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 34,
        "file": "Calcium-Ion.svg",
        "name": "Calcium Ion",
        "nameSlug": "calcium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 35,
        "file": "Calcium.svg",
        "name": "Calcium",
        "nameSlug": "calcium",
        "fileExtension": "svg"
    },
    {
        "ID": 36,
        "file": "Carbide.svg",
        "name": "Carbide",
        "nameSlug": "carbide",
        "fileExtension": "svg"
    },
    {
        "ID": 37,
        "file": "Carbonate.svg",
        "name": "Carbonate",
        "nameSlug": "carbonate",
        "fileExtension": "svg"
    },
    {
        "ID": 38,
        "file": "Carbon-Dioxide.svg",
        "name": "Carbon Dioxide",
        "nameSlug": "carbon-dioxide",
        "fileExtension": "svg"
    },
    {
        "ID": 39,
        "file": "Carbonic-Acid.svg",
        "name": "Carbonic Acid",
        "nameSlug": "carbonic-acid",
        "fileExtension": "svg"
    },
    {
        "ID": 40,
        "file": "Carbon-Ion.svg",
        "name": "Carbon Ion",
        "nameSlug": "carbon-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 41,
        "file": "Carbon-Monoxide.svg",
        "name": "Carbon Monoxide",
        "nameSlug": "carbon-monoxide",
        "fileExtension": "svg"
    },
    {
        "ID": 42,
        "file": "Carbon.svg",
        "name": "Carbon",
        "nameSlug": "carbon",
        "fileExtension": "svg"
    },
    {
        "ID": 43,
        "file": "Catalyst.svg",
        "name": "Catalyst",
        "nameSlug": "catalyst",
        "fileExtension": "svg"
    },
    {
        "ID": 44,
        "file": "Cesium-137.svg",
        "name": "Cesium 137",
        "nameSlug": "cesium-137",
        "fileExtension": "svg"
    },
    {
        "ID": 45,
        "file": "Chloride-Dots.svg",
        "name": "Chloride Dots",
        "nameSlug": "chloride-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 46,
        "file": "Chloride.svg",
        "name": "Chloride",
        "nameSlug": "chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 47,
        "file": "Chlorine-45.svg",
        "name": "Chlorine 45",
        "nameSlug": "chlorine-45",
        "fileExtension": "svg"
    },
    {
        "ID": 48,
        "file": "Chlorine-Atom.svg",
        "name": "Chlorine Atom",
        "nameSlug": "chlorine-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 49,
        "file": "Chlorine-Ion.svg",
        "name": "Chlorine Ion",
        "nameSlug": "chlorine-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 50,
        "file": "Chlorine-Nitrate.svg",
        "name": "Chlorine Nitrate",
        "nameSlug": "chlorine-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 51,
        "file": "Chlorine.svg",
        "name": "Chlorine",
        "nameSlug": "chlorine",
        "fileExtension": "svg"
    },
    {
        "ID": 52,
        "file": "Chromium-III-Oxide.svg",
        "name": "Chromium III Oxide",
        "nameSlug": "chromium-iii-oxide",
        "fileExtension": "svg"
    },
    {
        "ID": 53,
        "file": "Chromium-III.svg",
        "name": "Chromium III",
        "nameSlug": "chromium-iii",
        "fileExtension": "svg"
    },
    {
        "ID": 54,
        "file": "Chromium-II.svg",
        "name": "Chromium II",
        "nameSlug": "chromium-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 55,
        "file": "Chromium-Ion.svg",
        "name": "Chromium Ion",
        "nameSlug": "chromium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 56,
        "file": "Chromium.svg",
        "name": "Chromium",
        "nameSlug": "chromium",
        "fileExtension": "svg"
    },
    {
        "ID": 57,
        "file": "Copper-III.svg",
        "name": "Copper III",
        "nameSlug": "copper-iii",
        "fileExtension": "svg"
    },
    {
        "ID": 58,
        "file": "Copper-II-Nitrate.svg",
        "name": "Copper II Nitrate",
        "nameSlug": "copper-ii-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 59,
        "file": "Copper-II-Sulfate.svg",
        "name": "Copper II Sulfate",
        "nameSlug": "copper-ii-sulfate",
        "fileExtension": "svg"
    },
    {
        "ID": 60,
        "file": "Copper-II.svg",
        "name": "Copper II",
        "nameSlug": "copper-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 61,
        "file": "Copper.svg",
        "name": "Copper",
        "nameSlug": "copper",
        "fileExtension": "svg"
    },
    {
        "ID": 62,
        "file": "Cyanide-Dots.svg",
        "name": "Cyanide Dots",
        "nameSlug": "cyanide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 63,
        "file": "Cyanide.svg",
        "name": "Cyanide",
        "nameSlug": "cyanide",
        "fileExtension": "svg"
    },
    {
        "ID": 64,
        "file": "DDCF351F.png",
        "name": "DDCF351F",
        "nameSlug": "ddcf351f",
        "fileExtension": "png"
    },
    {
        "ID": 65,
        "file": "Deuterium.svg",
        "name": "Deuterium",
        "nameSlug": "deuterium",
        "fileExtension": "svg"
    },
    {
        "ID": 66,
        "file": "Dinitrogen-Tetroxide.svg",
        "name": "Dinitrogen Tetroxide",
        "nameSlug": "dinitrogen-tetroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 67,
        "file": "Energy.svg",
        "name": "Energy",
        "nameSlug": "energy",
        "fileExtension": "svg"
    },
    {
        "ID": 68,
        "file": "Ethanol.svg",
        "name": "Ethanol",
        "nameSlug": "ethanol",
        "fileExtension": "svg"
    },
    {
        "ID": 69,
        "file": "Ethene.svg",
        "name": "Ethene",
        "nameSlug": "ethene",
        "fileExtension": "svg"
    },
    {
        "ID": 70,
        "file": "filenames.txt",
        "name": "filenames.txt",
        "nameSlug": "filenames.txt",
        "fileExtension": "txt"
    },
    {
        "ID": 71,
        "file": "Fluoride.svg",
        "name": "Fluoride",
        "nameSlug": "fluoride",
        "fileExtension": "svg"
    },
    {
        "ID": 72,
        "file": "Fluorine-Atom.svg",
        "name": "Fluorine Atom",
        "nameSlug": "fluorine-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 73,
        "file": "Fluorine-Ion.svg",
        "name": "Fluorine Ion",
        "nameSlug": "fluorine-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 74,
        "file": "Fluorine.svg",
        "name": "Fluorine",
        "nameSlug": "fluorine",
        "fileExtension": "svg"
    },
    {
        "ID": 75,
        "file": "Generic.svg",
        "name": "Generic",
        "nameSlug": "generic",
        "fileExtension": "svg"
    },
    {
        "ID": 76,
        "file": "Glycerol.svg",
        "name": "Glycerol",
        "nameSlug": "glycerol",
        "fileExtension": "svg"
    },
    {
        "ID": 77,
        "file": "Gold-III.svg",
        "name": "Gold III",
        "nameSlug": "gold-iii",
        "fileExtension": "svg"
    },
    {
        "ID": 78,
        "file": "Gold-I.svg",
        "name": "Gold I",
        "nameSlug": "gold-i",
        "fileExtension": "svg"
    },
    {
        "ID": 79,
        "file": "Gold.svg",
        "name": "Gold",
        "nameSlug": "gold",
        "fileExtension": "svg"
    },
    {
        "ID": 80,
        "file": "Green.svg",
        "name": "Green",
        "nameSlug": "green",
        "fileExtension": "svg"
    },
    {
        "ID": 81,
        "file": "Helium-Ion.svg",
        "name": "Helium Ion",
        "nameSlug": "helium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 82,
        "file": "Helium.svg",
        "name": "Helium",
        "nameSlug": "helium",
        "fileExtension": "svg"
    },
    {
        "ID": 83,
        "file": "Hydrochloric-Acid-Dots.svg",
        "name": "Hydrochloric Acid Dots",
        "nameSlug": "hydrochloric-acid-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 84,
        "file": "Hydrochloric-Acid.svg",
        "name": "Hydrochloric Acid",
        "nameSlug": "hydrochloric-acid",
        "fileExtension": "svg"
    },
    {
        "ID": 85,
        "file": "Hydrogen-Atom.svg",
        "name": "Hydrogen Atom",
        "nameSlug": "hydrogen-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 86,
        "file": "Hydrogen-Bromide-Dots.svg",
        "name": "Hydrogen Bromide Dots",
        "nameSlug": "hydrogen-bromide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 87,
        "file": "Hydrogen-Bromide.svg",
        "name": "Hydrogen Bromide",
        "nameSlug": "hydrogen-bromide",
        "fileExtension": "svg"
    },
    {
        "ID": 88,
        "file": "Hydrogen-Cyanide-Dots.svg",
        "name": "Hydrogen Cyanide Dots",
        "nameSlug": "hydrogen-cyanide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 89,
        "file": "Hydrogen-Cyanide.svg",
        "name": "Hydrogen Cyanide",
        "nameSlug": "hydrogen-cyanide",
        "fileExtension": "svg"
    },
    {
        "ID": 90,
        "file": "Hydrogen-Fluoride.svg",
        "name": "Hydrogen Fluoride",
        "nameSlug": "hydrogen-fluoride",
        "fileExtension": "svg"
    },
    {
        "ID": 91,
        "file": "Hydrogen-Iodide.svg",
        "name": "Hydrogen Iodide",
        "nameSlug": "hydrogen-iodide",
        "fileExtension": "svg"
    },
    {
        "ID": 92,
        "file": "Hydrogen-Ion.svg",
        "name": "Hydrogen Ion",
        "nameSlug": "hydrogen-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 93,
        "file": "Hydrogen-Peroxide.svg",
        "name": "Hydrogen Peroxide",
        "nameSlug": "hydrogen-peroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 94,
        "file": "Hydrogen-Sulfide.svg",
        "name": "Hydrogen Sulfide",
        "nameSlug": "hydrogen-sulfide",
        "fileExtension": "svg"
    },
    {
        "ID": 95,
        "file": "Hydrogen.svg",
        "name": "Hydrogen",
        "nameSlug": "hydrogen",
        "fileExtension": "svg"
    },
    {
        "ID": 96,
        "file": "Hydronium.svg",
        "name": "Hydronium",
        "nameSlug": "hydronium",
        "fileExtension": "svg"
    },
    {
        "ID": 97,
        "file": "Hydroxide-Dots.svg",
        "name": "Hydroxide Dots",
        "nameSlug": "hydroxide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 98,
        "file": "Hydroxide.svg",
        "name": "Hydroxide",
        "nameSlug": "hydroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 99,
        "file": "Inert.svg",
        "name": "Inert",
        "nameSlug": "inert",
        "fileExtension": "svg"
    },
    {
        "ID": 100,
        "file": "Inhibitor.svg",
        "name": "Inhibitor",
        "nameSlug": "inhibitor",
        "fileExtension": "svg"
    },
    {
        "ID": 101,
        "file": "Iodine-Atom.svg",
        "name": "Iodine Atom",
        "nameSlug": "iodine-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 102,
        "file": "Iodine-Ion.svg",
        "name": "Iodine Ion",
        "nameSlug": "iodine-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 103,
        "file": "Iodine.svg",
        "name": "Iodine",
        "nameSlug": "iodine",
        "fileExtension": "svg"
    },
    {
        "ID": 104,
        "file": "Iron-III.svg",
        "name": "Iron III",
        "nameSlug": "iron-iii",
        "fileExtension": "svg"
    },
    {
        "ID": 105,
        "file": "Iron-II-Sulfate.svg",
        "name": "Iron II Sulfate",
        "nameSlug": "iron-ii-sulfate",
        "fileExtension": "svg"
    },
    {
        "ID": 106,
        "file": "Iron-II.svg",
        "name": "Iron II",
        "nameSlug": "iron-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 107,
        "file": "Iron-I.svg",
        "name": "Iron I",
        "nameSlug": "iron-i",
        "fileExtension": "svg"
    },
    {
        "ID": 108,
        "file": "Iron.svg",
        "name": "Iron",
        "nameSlug": "iron",
        "fileExtension": "svg"
    },
    {
        "ID": 109,
        "file": "Krypton-50.svg",
        "name": "Krypton 50",
        "nameSlug": "krypton-50",
        "fileExtension": "svg"
    },
    {
        "ID": 110,
        "file": "Krypton-89.svg",
        "name": "Krypton 89",
        "nameSlug": "krypton-89",
        "fileExtension": "svg"
    },
    {
        "ID": 111,
        "file": "Lead-Iodide.svg",
        "name": "Lead Iodide",
        "nameSlug": "lead-iodide",
        "fileExtension": "svg"
    },
    {
        "ID": 112,
        "file": "Lead-Ion.svg",
        "name": "Lead Ion",
        "nameSlug": "lead-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 113,
        "file": "Lead.svg",
        "name": "Lead",
        "nameSlug": "lead",
        "fileExtension": "svg"
    },
    {
        "ID": 114,
        "file": "Lithium-Chloride.svg",
        "name": "Lithium Chloride",
        "nameSlug": "lithium-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 115,
        "file": "Lithium-Hydroxide.svg",
        "name": "Lithium Hydroxide",
        "nameSlug": "lithium-hydroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 116,
        "file": "Lithium-Ion.svg",
        "name": "Lithium Ion",
        "nameSlug": "lithium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 117,
        "file": "Lithium-I.svg",
        "name": "Lithium I",
        "nameSlug": "lithium-i",
        "fileExtension": "svg"
    },
    {
        "ID": 118,
        "file": "Lithium-Nitrate.svg",
        "name": "Lithium Nitrate",
        "nameSlug": "lithium-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 119,
        "file": "Lithium-Sulfide.svg",
        "name": "Lithium Sulfide",
        "nameSlug": "lithium-sulfide",
        "fileExtension": "svg"
    },
    {
        "ID": 120,
        "file": "Lithium.svg",
        "name": "Lithium",
        "nameSlug": "lithium",
        "fileExtension": "svg"
    },
    {
        "ID": 121,
        "file": "Magnesium-Ion.svg",
        "name": "Magnesium Ion",
        "nameSlug": "magnesium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 122,
        "file": "Magnesium-Sulfate.svg",
        "name": "Magnesium Sulfate",
        "nameSlug": "magnesium-sulfate",
        "fileExtension": "svg"
    },
    {
        "ID": 123,
        "file": "Magnesium.svg",
        "name": "Magnesium",
        "nameSlug": "magnesium",
        "fileExtension": "svg"
    },
    {
        "ID": 124,
        "file": "Manganese-Dioxide.svg",
        "name": "Manganese Dioxide",
        "nameSlug": "manganese-dioxide",
        "fileExtension": "svg"
    },
    {
        "ID": 125,
        "file": "Manganese-II-Chloride.svg",
        "name": "Manganese II Chloride",
        "nameSlug": "manganese-ii-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 126,
        "file": "Manganese-II.svg",
        "name": "Manganese II",
        "nameSlug": "manganese-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 127,
        "file": "Manganese-IV.svg",
        "name": "Manganese IV",
        "nameSlug": "manganese-iv",
        "fileExtension": "svg"
    },
    {
        "ID": 128,
        "file": "Manganese.svg",
        "name": "Manganese",
        "nameSlug": "manganese",
        "fileExtension": "svg"
    },
    {
        "ID": 129,
        "file": "Mercury-II.svg",
        "name": "Mercury II",
        "nameSlug": "mercury-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 130,
        "file": "Mercury-I.svg",
        "name": "Mercury I",
        "nameSlug": "mercury-i",
        "fileExtension": "svg"
    },
    {
        "ID": 131,
        "file": "Mercury.svg",
        "name": "Mercury",
        "nameSlug": "mercury",
        "fileExtension": "svg"
    },
    {
        "ID": 132,
        "file": "Methane.svg",
        "name": "Methane",
        "nameSlug": "methane",
        "fileExtension": "svg"
    },
    {
        "ID": 133,
        "file": "Methylamine.svg",
        "name": "Methylamine",
        "nameSlug": "methylamine",
        "fileExtension": "svg"
    },
    {
        "ID": 134,
        "file": "Methylammonium.svg",
        "name": "Methylammonium",
        "nameSlug": "methylammonium",
        "fileExtension": "svg"
    },
    {
        "ID": 135,
        "file": "Neutron-small.svg",
        "name": "Neutron small",
        "nameSlug": "neutron-small",
        "fileExtension": "svg"
    },
    {
        "ID": 136,
        "file": "Neutron.svg",
        "name": "Neutron",
        "nameSlug": "neutron",
        "fileExtension": "svg"
    },
    {
        "ID": 137,
        "file": "Nitrate.svg",
        "name": "Nitrate",
        "nameSlug": "nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 138,
        "file": "Nitric-Acid.svg",
        "name": "Nitric Acid",
        "nameSlug": "nitric-acid",
        "fileExtension": "svg"
    },
    {
        "ID": 139,
        "file": "Nitric-Oxide.svg",
        "name": "Nitric Oxide",
        "nameSlug": "nitric-oxide",
        "fileExtension": "svg"
    },
    {
        "ID": 140,
        "file": "Nitrogen-Atom.svg",
        "name": "Nitrogen Atom",
        "nameSlug": "nitrogen-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 141,
        "file": "Nitrogen-Dioxide.svg",
        "name": "Nitrogen Dioxide",
        "nameSlug": "nitrogen-dioxide",
        "fileExtension": "svg"
    },
    {
        "ID": 142,
        "file": "Nitrogen-Ion.svg",
        "name": "Nitrogen Ion",
        "nameSlug": "nitrogen-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 143,
        "file": "Nitrogen.svg",
        "name": "Nitrogen",
        "nameSlug": "nitrogen",
        "fileExtension": "svg"
    },
    {
        "ID": 144,
        "file": "Nitrosyl-Chloride.svg",
        "name": "Nitrosyl Chloride",
        "nameSlug": "nitrosyl-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 145,
        "file": "Nitryl-Chloride.svg",
        "name": "Nitryl Chloride",
        "nameSlug": "nitryl-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 146,
        "file": "Nuclear-01.svg",
        "name": "Nuclear 01",
        "nameSlug": "nuclear-01",
        "fileExtension": "svg"
    },
    {
        "ID": 147,
        "file": "Nuclear-02.svg",
        "name": "Nuclear 02",
        "nameSlug": "nuclear-02",
        "fileExtension": "svg"
    },
    {
        "ID": 148,
        "file": "Nuclear-03.svg",
        "name": "Nuclear 03",
        "nameSlug": "nuclear-03",
        "fileExtension": "svg"
    },
    {
        "ID": 149,
        "file": "Nuclear-04.svg",
        "name": "Nuclear 04",
        "nameSlug": "nuclear-04",
        "fileExtension": "svg"
    },
    {
        "ID": 150,
        "file": "Nuclear-05.svg",
        "name": "Nuclear 05",
        "nameSlug": "nuclear-05",
        "fileExtension": "svg"
    },
    {
        "ID": 151,
        "file": "Nuclear-06.svg",
        "name": "Nuclear 06",
        "nameSlug": "nuclear-06",
        "fileExtension": "svg"
    },
    {
        "ID": 152,
        "file": "Nuclear-07.svg",
        "name": "Nuclear 07",
        "nameSlug": "nuclear-07",
        "fileExtension": "svg"
    },
    {
        "ID": 153,
        "file": "Nuclear-p19-01.svg",
        "name": "Nuclear p19 01",
        "nameSlug": "nuclear-p19-01",
        "fileExtension": "svg"
    },
    {
        "ID": 154,
        "file": "Nuclear-p22-01.svg",
        "name": "Nuclear p22 01",
        "nameSlug": "nuclear-p22-01",
        "fileExtension": "svg"
    },
    {
        "ID": 155,
        "file": "Nuclear-p22-05.svg",
        "name": "Nuclear p22 05",
        "nameSlug": "nuclear-p22-05",
        "fileExtension": "svg"
    },
    {
        "ID": 156,
        "file": "Nuclear-p39-01.svg",
        "name": "Nuclear p39 01",
        "nameSlug": "nuclear-p39-01",
        "fileExtension": "svg"
    },
    {
        "ID": 157,
        "file": "Nuclear-p39-02.svg",
        "name": "Nuclear p39 02",
        "nameSlug": "nuclear-p39-02",
        "fileExtension": "svg"
    },
    {
        "ID": 158,
        "file": "Octane.svg",
        "name": "Octane",
        "nameSlug": "octane",
        "fileExtension": "svg"
    },
    {
        "ID": 159,
        "file": "Oxygen-Atom.svg",
        "name": "Oxygen Atom",
        "nameSlug": "oxygen-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 160,
        "file": "Oxygen-Ion.svg",
        "name": "Oxygen Ion",
        "nameSlug": "oxygen-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 161,
        "file": "Oxygen.svg",
        "name": "Oxygen",
        "nameSlug": "oxygen",
        "fileExtension": "svg"
    },
    {
        "ID": 162,
        "file": "Ozone.svg",
        "name": "Ozone",
        "nameSlug": "ozone",
        "fileExtension": "svg"
    },
    {
        "ID": 163,
        "file": "Particle-BetaNegative.svg",
        "name": "Particle BetaNegative",
        "nameSlug": "particle-betanegative",
        "fileExtension": "svg"
    },
    {
        "ID": 164,
        "file": "Pentane.svg",
        "name": "Pentane",
        "nameSlug": "pentane",
        "fileExtension": "svg"
    },
    {
        "ID": 165,
        "file": "Phenylpthalein.svg",
        "name": "Phenylpthalein",
        "nameSlug": "phenylpthalein",
        "fileExtension": "svg"
    },
    {
        "ID": 166,
        "file": "Phosphorus-Ion.svg",
        "name": "Phosphorus Ion",
        "nameSlug": "phosphorus-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 167,
        "file": "Phosphorus-Pentachloride.svg",
        "name": "Phosphorus Pentachloride",
        "nameSlug": "phosphorus-pentachloride",
        "fileExtension": "svg"
    },
    {
        "ID": 168,
        "file": "Phosphorus.svg",
        "name": "Phosphorus",
        "nameSlug": "phosphorus",
        "fileExtension": "svg"
    },
    {
        "ID": 169,
        "file": "Phosphorus-Trichloride.svg",
        "name": "Phosphorus Trichloride",
        "nameSlug": "phosphorus-trichloride",
        "fileExtension": "svg"
    },
    {
        "ID": 170,
        "file": "Potassium-85.svg",
        "name": "Potassium 85",
        "nameSlug": "potassium-85",
        "fileExtension": "svg"
    },
    {
        "ID": 171,
        "file": "Potassium-Bromide.svg",
        "name": "Potassium Bromide",
        "nameSlug": "potassium-bromide",
        "fileExtension": "svg"
    },
    {
        "ID": 172,
        "file": "Potassium-Chloride.svg",
        "name": "Potassium Chloride",
        "nameSlug": "potassium-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 173,
        "file": "Potassium-Ion.svg",
        "name": "Potassium Ion",
        "nameSlug": "potassium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 174,
        "file": "Potassium-Nitrate.svg",
        "name": "Potassium Nitrate",
        "nameSlug": "potassium-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 175,
        "file": "Potassium.svg",
        "name": "Potassium",
        "nameSlug": "potassium",
        "fileExtension": "svg"
    },
    {
        "ID": 176,
        "file": "Propane.svg",
        "name": "Propane",
        "nameSlug": "propane",
        "fileExtension": "svg"
    },
    {
        "ID": 177,
        "file": "Proton-small.svg",
        "name": "Proton small",
        "nameSlug": "proton-small",
        "fileExtension": "svg"
    },
    {
        "ID": 178,
        "file": "Proton.svg",
        "name": "Proton",
        "nameSlug": "proton",
        "fileExtension": "svg"
    },
    {
        "ID": 179,
        "file": "Quark-down-noShadow.svg",
        "name": "Quark down noShadow",
        "nameSlug": "quark-down-noshadow",
        "fileExtension": "svg"
    },
    {
        "ID": 180,
        "file": "Quark-down.svg",
        "name": "Quark down",
        "nameSlug": "quark-down",
        "fileExtension": "svg"
    },
    {
        "ID": 181,
        "file": "Quark-up-noShadow.svg",
        "name": "Quark up noShadow",
        "nameSlug": "quark-up-noshadow",
        "fileExtension": "svg"
    },
    {
        "ID": 182,
        "file": "Quark-up.svg",
        "name": "Quark up",
        "nameSlug": "quark-up",
        "fileExtension": "svg"
    },
    {
        "ID": 183,
        "file": "Red.svg",
        "name": "Red",
        "nameSlug": "red",
        "fileExtension": "svg"
    },
    {
        "ID": 184,
        "file": "Rubidium-96.svg",
        "name": "Rubidium 96",
        "nameSlug": "rubidium-96",
        "fileExtension": "svg"
    },
    {
        "ID": 185,
        "file": "Shape-1.svg",
        "name": "Shape 1",
        "nameSlug": "shape-1",
        "fileExtension": "svg"
    },
    {
        "ID": 186,
        "file": "Shape-2.svg",
        "name": "Shape 2",
        "nameSlug": "shape-2",
        "fileExtension": "svg"
    },
    {
        "ID": 187,
        "file": "Silicon-Dioxide.svg",
        "name": "Silicon Dioxide",
        "nameSlug": "silicon-dioxide",
        "fileExtension": "svg"
    },
    {
        "ID": 188,
        "file": "Silicon-Ion.svg",
        "name": "Silicon Ion",
        "nameSlug": "silicon-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 189,
        "file": "Silicon-IV.svg",
        "name": "Silicon IV",
        "nameSlug": "silicon-iv",
        "fileExtension": "svg"
    },
    {
        "ID": 190,
        "file": "Silicon.svg",
        "name": "Silicon",
        "nameSlug": "silicon",
        "fileExtension": "svg"
    },
    {
        "ID": 191,
        "file": "Silver-Bromide.svg",
        "name": "Silver Bromide",
        "nameSlug": "silver-bromide",
        "fileExtension": "svg"
    },
    {
        "ID": 192,
        "file": "Silver-Carbonate.svg",
        "name": "Silver Carbonate",
        "nameSlug": "silver-carbonate",
        "fileExtension": "svg"
    },
    {
        "ID": 193,
        "file": "Silver-Chloride.svg",
        "name": "Silver Chloride",
        "nameSlug": "silver-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 194,
        "file": "Silver-Hydroxide.svg",
        "name": "Silver Hydroxide",
        "nameSlug": "silver-hydroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 195,
        "file": "Silver-Ion.svg",
        "name": "Silver Ion",
        "nameSlug": "silver-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 196,
        "file": "Silver-Nitrate.svg",
        "name": "Silver Nitrate",
        "nameSlug": "silver-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 197,
        "file": "Silver.svg",
        "name": "Silver",
        "nameSlug": "silver",
        "fileExtension": "svg"
    },
    {
        "ID": 198,
        "file": "Sodium-Acetate.svg",
        "name": "Sodium Acetate",
        "nameSlug": "sodium-acetate",
        "fileExtension": "svg"
    },
    {
        "ID": 199,
        "file": "Sodium-Bicarbonate.svg",
        "name": "Sodium Bicarbonate",
        "nameSlug": "sodium-bicarbonate",
        "fileExtension": "svg"
    },
    {
        "ID": 200,
        "file": "Sodium-Carbonate.svg",
        "name": "Sodium Carbonate",
        "nameSlug": "sodium-carbonate",
        "fileExtension": "svg"
    },
    {
        "ID": 201,
        "file": "Sodium-Chloride-Dots.svg",
        "name": "Sodium Chloride Dots",
        "nameSlug": "sodium-chloride-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 202,
        "file": "Sodium-Chloride.svg",
        "name": "Sodium Chloride",
        "nameSlug": "sodium-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 203,
        "file": "Sodium-Hydroxide-Dots.svg",
        "name": "Sodium Hydroxide Dots",
        "nameSlug": "sodium-hydroxide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 204,
        "file": "Sodium-Hydroxide.svg",
        "name": "Sodium Hydroxide",
        "nameSlug": "sodium-hydroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 205,
        "file": "Sodium-Ion.svg",
        "name": "Sodium Ion",
        "nameSlug": "sodium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 206,
        "file": "Sodium-I.svg",
        "name": "Sodium I",
        "nameSlug": "sodium-i",
        "fileExtension": "svg"
    },
    {
        "ID": 207,
        "file": "Sodium-Nitrate.svg",
        "name": "Sodium Nitrate",
        "nameSlug": "sodium-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 208,
        "file": "Sodium-Nitrite.svg",
        "name": "Sodium Nitrite",
        "nameSlug": "sodium-nitrite",
        "fileExtension": "svg"
    },
    {
        "ID": 209,
        "file": "Sodium.svg",
        "name": "Sodium",
        "nameSlug": "sodium",
        "fileExtension": "svg"
    },
    {
        "ID": 210,
        "file": "Spark.svg",
        "name": "Spark",
        "nameSlug": "spark",
        "fileExtension": "svg"
    },
    {
        "ID": 211,
        "file": "Sulfate.svg",
        "name": "Sulfate",
        "nameSlug": "sulfate",
        "fileExtension": "svg"
    },
    {
        "ID": 212,
        "file": "Sulfide.svg",
        "name": "Sulfide",
        "nameSlug": "sulfide",
        "fileExtension": "svg"
    },
    {
        "ID": 213,
        "file": "Sulfur-Dioxide.svg",
        "name": "Sulfur Dioxide",
        "nameSlug": "sulfur-dioxide",
        "fileExtension": "svg"
    },
    {
        "ID": 214,
        "file": "Sulfur-Ion.svg",
        "name": "Sulfur Ion",
        "nameSlug": "sulfur-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 215,
        "file": "Sulfur.svg",
        "name": "Sulfur",
        "nameSlug": "sulfur",
        "fileExtension": "svg"
    },
    {
        "ID": 216,
        "file": "Sulfur-Tetrafluoride.svg",
        "name": "Sulfur Tetrafluoride",
        "nameSlug": "sulfur-tetrafluoride",
        "fileExtension": "svg"
    },
    {
        "ID": 217,
        "file": "Thiocyanate.svg",
        "name": "Thiocyanate",
        "nameSlug": "thiocyanate",
        "fileExtension": "svg"
    },
    {
        "ID": 218,
        "file": "Thiocyanic-Acid.svg",
        "name": "Thiocyanic Acid",
        "nameSlug": "thiocyanic-acid",
        "fileExtension": "svg"
    },
    {
        "ID": 219,
        "file": "Tin-II.svg",
        "name": "Tin II",
        "nameSlug": "tin-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 220,
        "file": "Tin-IV.svg",
        "name": "Tin IV",
        "nameSlug": "tin-iv",
        "fileExtension": "svg"
    },
    {
        "ID": 221,
        "file": "Tin.svg",
        "name": "Tin",
        "nameSlug": "tin",
        "fileExtension": "svg"
    },
    {
        "ID": 222,
        "file": "trichlorofluoromethane.svg",
        "name": "Trichlorofluoromethane",
        "nameSlug": "trichlorofluoromethane",
        "fileExtension": "svg"
    },
    {
        "ID": 223,
        "file": "Uranium-235.svg",
        "name": "Uranium-235",
        "nameSlug": "uranium-235",
        "fileExtension": "svg"
    },
    {
        "ID": 224,
        "file": "Water-Dots.svg",
        "name": "Water (Lewis)",
        "nameSlug": "water-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 225,
        "file": "Water.svg",
        "name": "Water",
        "nameSlug": "water",
        "fileExtension": "svg"
    },
    {
        "ID": 226,
        "file": "Zinc-Ion.svg",
        "name": "Zinc Ion",
        "nameSlug": "zinc-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 227,
        "file": "Zinc.svg",
        "name": "Zinc",
        "nameSlug": "zinc",
        "fileExtension": "svg"
    }
]