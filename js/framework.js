var gl, gl_draw_buffers;
var width, height;

(function() {
    'use strict';

    var canvas, renderer, scene, camera, controls, stats;
    var models = [];

    var cameraMat = new THREE.Matrix4();

    var render = function() {
        if (typeof R.gridInfo !== 'undefined'){
            camera.updateMatrixWorld();
            camera.matrixWorldInverse.getInverse(camera.matrixWorld);
            cameraMat.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);

            R.particleRender({
                models: models,
                cameraMat: cameraMat,
                camera: camera,
                cameraPos: camera.getWorldPosition()
            });
        }
    };

    R.update = function() {
        controls.update();
        stats.end();
        stats.begin();
        if (R.toReset) {
            R.time = 0;
            R.toReset = false;
            R.setupBuffers();
        }
        render();
        if (!aborted) {
            requestAnimationFrame(R.update);
        }
    };

    var resize = function() {
        camera.aspect = width / height;
        camera.updateProjectionMatrix();
        renderer.setSize(width, height);
        render();
    };

    var initExtensions = function() {
        var extensions = gl.getSupportedExtensions();
        console.log(extensions);

        var reqd = [
            'OES_texture_float',
            'OES_texture_float_linear',
            'WEBGL_depth_texture',
            'WEBGL_draw_buffers',
            'EXT_frag_depth'
        ];
        for (var i = 0; i < reqd.length; i++) {
            var e = reqd[i];
            if (extensions.indexOf(e) < 0) {
                abort('unable to load extension: ' + e);
            }
        }

        gl.getExtension('OES_texture_float');
        gl.getExtension('OES_texture_float_linear');
        gl.getExtension('WEBGL_depth_texture');
        gl.getExtension('EXT_frag_depth');

        gl_draw_buffers = gl.getExtension('WEBGL_draw_buffers');
        var maxdb = gl.getParameter(gl_draw_buffers.MAX_DRAW_BUFFERS_WEBGL);
        console.log('MAX_DRAW_BUFFERS_WEBGL: ' + maxdb);
    };

    var init = function() {
        var debugMode = false;

        canvas = document.getElementById('canvas');
        renderer = new THREE.WebGLRenderer({
            canvas: canvas,
            preserveDrawingBuffer: debugMode
        });
        gl = renderer.context;

        if (debugMode) {
            $('#debugmodewarning').css('display', 'block');
            var throwOnGLError = function(err, funcName, args) {
                abort(WebGLDebugUtils.glEnumToString(err) +
                    " was caused by call to: " + funcName);
            };
            gl = WebGLDebugUtils.makeDebugContext(gl, throwOnGLError);
        }

        initExtensions();

        stats = new Stats();
        stats.setMode(1); // 0: fps, 1: ms, 2: mb
        stats.domElement.style.position = 'absolute';
        stats.domElement.style.left = '0px';
        stats.domElement.style.top = '0px';
        document.body.appendChild(stats.domElement);

        scene = new THREE.Scene();

        width = canvas.width;
        height = canvas.height;
        R.fovy = 45;
        camera = new THREE.PerspectiveCamera(
            R.fovy,             // Field of view
            width / height, // Aspect ratio
            0.1,            // Near plane
            100             // Far plane
        );
        //camera.position.set(-3, 3, -3);
        camera.position.set(0, 3, 3.5);
        //camera.position.set(0, 4, 1.5);
        R.nearPlaneHeight = height / (2*Math.tan(0.5* R.fovy*Math.PI/180.0));
        //console.log(nearPlaneHeight);

        controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.enableZoom = true;
        //controls.target.set(0.5, 0.5, 0.5);
        //controls.target.set(0.0, -0.1, 0.0);
        controls.target.set(0.0, 0.4, 0.0);
        controls.rotateSpeed = 0.3;
        controls.zoomSpeed = 1.0;
        controls.panSpeed = 2.0;

        R.model = models[0];
	    R.initBuffers();
	    R.loadAllShaderPrograms();
	    requestAnimationFrame(R.update);
        resize();

        gl.clearColor(0.5, 0.5, 0.5, 0.5);
        gl.clearDepth(1.0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);        

    };

    window.handle_load.push(init);
})();
