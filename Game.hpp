#pragma once

#include "GL.hpp"

#include <SDL.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include <vector>
#include <random>

// The 'Game' struct holds all of the game-relevant state,
// and is called by the main loop.

struct Game {
	//Game creates OpenGL resources (i.e. vertex buffer objects) in its
	//constructor and frees them in its destructor.
	Game();
	~Game();

	//handle_event is called when new mouse or keyboard events are received:
	// (note that this might be many times per frame or never)
	//The function should return 'true' if it handled the event.
	bool handle_event(SDL_Event const &evt, glm::uvec2 window_size);

	//update is called at the start of a new frame, after events are handled:
	void update(float elapsed);

	//draw is called after update:
	void draw(glm::uvec2 drawable_size);

  //game logic functions
  glm::uvec2 move_by(glm::uvec2 pos, int fx, int fy);

  void move(int x, int y);

  void attack();

  void next_turn();

  bool check_goblin_collisions();

  void move_goblins();

  bool space_free(glm::uint x, glm::uint y);

  void restart();

  bool check_win();

  void generate_level();

  glm::ivec2 get_goblin_move(glm::uvec2 g_pos, glm::uvec2 p_pos);

  bool is_valid_space(glm::uint px, glm::uint py, glm::ivec2 facing);

  glm::ivec2 get_random_facing(unsigned int seed);

	//------- opengl resources -------

	//shader program that draws lit objects with vertex colors:
	struct {
		GLuint program = -1U; //program object

		//uniform locations:
		GLuint object_to_clip_mat4 = -1U;
		GLuint object_to_light_mat4x3 = -1U;
		GLuint normal_to_light_mat3 = -1U;
		GLuint sun_direction_vec3 = -1U;
		GLuint sun_color_vec3 = -1U;
		GLuint sky_direction_vec3 = -1U;
		GLuint sky_color_vec3 = -1U;

		//attribute locations:
		GLuint Position_vec4 = -1U;
		GLuint Normal_vec3 = -1U;
		GLuint Color_vec4 = -1U;
	} simple_shading;

	//mesh data, stored in a vertex buffer:
	GLuint meshes_vbo = -1U; //vertex buffer holding mesh data

	//The location of each mesh in the meshes vertex buffer:
	struct Mesh {
		GLint first = 0;
		GLsizei count = 0;
	};

  Mesh player_mesh;
  Mesh goblin_mesh;
  Mesh target_mesh;
  Mesh bg_mesh;

	GLuint meshes_for_simple_shading_vao = -1U; //vertex array object that describes how to connect the meshes_vbo to the simple_shading_program

	//------- game state -------

	glm::uvec2 board_size = glm::uvec2(6,6);
	std::vector< Mesh const * > board_meshes;
	std::vector< glm::quat > board_rotations;

	glm::uvec2 player_pos = glm::uvec2(0, 0);
  glm::ivec2 player_facing = glm::ivec2(1, 0);
  glm::uvec2 player_start_pos = glm::uvec2(0, 0);
  glm::ivec2 player_start_facing = glm::ivec2(1, 0);

  std::vector< glm::uvec2 > goblin_positions;
  std::vector< glm::uvec2 > goblin_start_positions;

  std::mt19937 mt;
};
