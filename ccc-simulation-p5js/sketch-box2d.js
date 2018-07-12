// Molecule variables
var moleculeArray = [];
var collectionArray = [];
var database;

// Boundaries
var boundary;
var piston;

// System variables
var systemKineticEnergy;

// Global settings

// Sets the relative directory to the folder containing particle sprites
var imageDir = "img/svg/";

// Scale factor used to scale distances and velocities in b2
var b2ScaleFactor = 3;

// Index to assign each particle a unique ID
var particleIndex = 0;

// Default force constant for intermolecular interactions
var defaultIMForceConstant = 1e2;

function preload() {
    // Database is read in from an external JSON
    database = loadJSON("assets/database.json");

    // Reads in simulation configuration file from external settings JSON
    settings = loadJSON("assets/sample-config.json");
}

function setup() {
    // Creates the drawing canvas
    createCanvas(800, 600);

    // Sets all x,y coordinates to the center of images
    imageMode(CENTER);

    // Adjust framerate to slow animation, useful for debugging
    // frameRate(5);

    // Initiates a new b2 world, scaling factor, gravity vector
    b2newWorld(settings.global.b2ScaleFactor || b2ScaleFactor, createVector(0, 9.8));

    // Creates the boundary for a closed reaction container
    boundary = new Boundary("CLOSED");

    // Creates the piston boundary for the reaction container
    piston = new Boundary("PISTON", createVector(width / 2, 0), createVector(width, 20))

    // Creates a collection of particles from settings file
    for (var i = 0; i < Object.keys(settings.collections).length; i++) {
        collectionArray[i] = new Collection(settings.collections[i].key, createVector(settings.collections[i].position[0], settings.collections[i].position[1]), settings.collections[i].separation, settings.collections[i].numberOfParticles, settings.collections[i].columns, settings.collections[i].options);
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
    for (var i = 0; i < moleculeArray.length; i++) {
        for (var j = 0; j < moleculeArray.length; j++) {
            if (moleculeArray[i].id == moleculeArray[j].id || moleculeArray[i].indexedBy(moleculeArray[j])) {
                continue;
            } else {
                // Stores that molecule_i and molecule_j are a known force pair
                // Prevents duplicate calculation of IM forces per each draw loop
                moleculeArray[i].indexes(moleculeArray[j]);

                // Apply the force from the vdW interactions (using a Coulombic potential function)
                // Currently the code assumes that when dissimilar particles interact, that the lesser of the two IM force constants should be used
                var imVector = moleculeArray[i].calculateIMForce(moleculeArray[j], Math.min(moleculeArray[i].getIMForceConstant(), moleculeArray[j].getIMForceConstant()));

                // Newton's 3rd Law
                moleculeArray[i].addForceToNetForce(imVector);
                moleculeArray[j].addForceToNetForce(imVector.mult(-1));

            }
        }
    }

    // Iterates over particles to update system KE & apply forces
    for (var i = 0; i < moleculeArray.length; i++) {
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

// Removed the setup logic for the initial simulation of a default set of particles
function defaultParticles(numberOfParticles) {
    // Calculates how many rows of particles to draw, currently assumes 5 columns of particles
    var rowsOfParticles = Math.floor(numberOfParticles / 5)

    // Increments the number of rows if there isn't exactly a multiple of 5 particles
    if (numberOfParticles % 5 > 0) {
        rowsOfParticles++;
    }

    // Creates molecule collection in a grid with random velocity vectors
    for (var i = 0; i < rowsOfParticles; i++) {
        for (var j = 0; j < 5; j++) {
            moleculeArray[5 * i + j] = new Particle(225, createVector(width / 10 + (j * (width / 5)), height / 10 + i * (height / 10)), p5.Vector.random2D().mult(random(1, 5)), 5 * i + j);
        }
    }
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

    constructor(key, position, velocity, id, options) {
        // This is the value in the "ID" field of the database variable
        this.databaseKey = key;

        // This gives each molecule a unique index to track it within a set of molecules
        this.id = id;

        // Particle ID is database key, used to pull name
        this.name = database[this.databaseKey - 1].name;

        // Sets the position of the particle
        this.position = createVector(position.x, position.y);

        // Sets the velocity of the particle
        this.velocity = createVector(velocity.x, velocity.y);

        // Defines the intermolecular force scale constant
        this.imForceConstant = options.imForceConstant || defaultIMForceConstant;

        // Gravity scalar
        this.gravityScale = options.gravityScale || 1;

        // Initializes the particle mass
        this.mass = database[this.databaseKey - 1].mass || 1;

        // When forces are applied, we need to indicate which particle force pairs have already been calculated
        // This avoids duplicate calculation of intermolecular forces
        this.forceIndices = [];

        // Net force vector on the particle
        this.netForce = createVector(0, 0);

        // Sets the initial KE of the particle, we start at rest before applying an initial velocity, therefore KE = 0
        this.kineticEnergy = 0;

        // Sets the image reference for the molecule using the ID number of the particle
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
            var fOverR = 12 * epsilonValue * (Math.pow((this.getRadius() / distanceBetweenParticles), 13) - Math.pow((this.getRadius() / distanceBetweenParticles), 7)) * Math.pow(distanceBetweenParticles, -1);

            forceToReturn = Math.min(fOverR * distanceBetweenParticles, forceCutoffValue);
        }
        // Negative forces are attractive by definition, so we multiply by -1 to ensure the force vector has the same direction as r
        return vectorBetweenParticles.normalize().mult(-forceToReturn);
    }

    // A simplified 1/r^2 force equation
    calculateIMForce(particle, scaleFactor) {
        var forcetoReturn, scaleFactorValue;
        var distanceBetweenParticles = this.distanceToParticle(particle);
        var vectorBetweenParticles = this.vectorToParticle(particle);

        if (typeof scaleFactor == "undefined") {
            scaleFactorValue = 1.0;
        } else {
            scaleFactorValue = scaleFactor;
        }

        // Imposes a efficiency constraint to limit the number of calculations per cycle
        if ((distanceBetweenParticles / this.getDiameter()) > 3) {
            forcetoReturn = 0;
        } else {
            forcetoReturn = scaleFactorValue * Math.pow((distanceBetweenParticles / this.getDiameter()), 2);
        }

        return vectorBetweenParticles.normalize().mult(forcetoReturn);
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

    getDiameter() {
        return (!(typeof this.body == 'undefined')) ? Math.max(this.body.wh(0).x, this.body.wh(0).y) : Math.max(this.imageObject.width, this.imageObject.height);
    }

    // Returns the net force vector (intermolecular forces)
    getNetForce() {
        return this.netForce;
    }

    // Returns the Intermolecular force constant
    getIMForceConstant() {
        return this.imForceConstant;
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
        this.body = new b2Body('box', true, createVector(this.position.x, this.position.y), createVector(this.imageObject.width, this.imageObject.height), 1.0, 0, 1.0);
        this.body.image(imageObject, 0);

        // Sets the gravity scalar per particle
        this.body.gravityScale = this.gravityScale;

        // Set particle velocity
        this.body.applyForce(createVector(this.velocity.x, this.velocity.y), 1e4);

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

// Collection of particles - used to create grids of particles
class Collection {
    constructor(key, position, separation, numberOfParticles, columns, options) {
        // This is the value in the "ID" field of the database variable
        this.databaseKey = key;

        // Sets the default position of the top left atom
        this.position = createVector(position.x, position.y);

        // Sets the default separation between particles in pixels
        this.separation = separation || 0;

        // Sets how many particles are in this collection
        this.numberOfParticles = numberOfParticles || 1;

        // Defines how many columns the particles will be drawn in
        this.columns = columns || 1;

        // Sets the force constant for the particles
        this.options = options;

        // Begins the drawing of particles
        this.initializeParticles();
    }

    initializeParticles() {
        // Calculates how many rows of particles to draw, currently assumes 5 columns of particles
        var rowsOfParticles = Math.floor(this.numberOfParticles / this.columns);

        // Increments the number of rows if there isn't exactly a multiple of 5 particles
        if (this.numberOfParticles % this.columns > 0) {
            rowsOfParticles++;
        }

        // Creates molecule collection in a grid with random velocity vectors
        for (var i = 0; i < rowsOfParticles; i++) {
            for (var j = 0; j < this.columns; j++) {
                
                // Adds the particle to the existing molecule array
                // This is just a quick implementation - more robust would keep particles organized in collection
                moleculeArray[particleIndex] = new Particle(this.databaseKey, createVector(this.position.x + (j * this.separation), this.position.y + i * this.separation), createVector(0,0), particleIndex, this.options);

                // Increments the global particle index so that unique values are used
                particleIndex++;
            }
        }
    }


}

