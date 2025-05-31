
#pragma once
#include <cstdint>
#include <SFML/System/Vector2.hpp>

using Vec2  = sf::Vector2f;
using IVec2 = sf::Vector2i;
#include <vector>
#include <array>


template<typename T>
struct Grid
{
	struct HitPoint
	{
		T* cell;
		float dist;

		HitPoint()
			: cell(nullptr)
			, dist(0.0f)
		{}
	};

	int32_t width, height;
	std::vector<T> data;

	Grid()
		: width(0)
		, height(0)
	{}

	Grid(int32_t width_, int32_t height_)
		: width(width_)
		, height(height_)
	{
		data.resize(width * height);
	}

	int32_t mod(int32_t dividend, int32_t divisor) const
	{
		return (dividend%divisor + divisor) % divisor;
	}

	template<typename Vec2Type>
	bool checkCoords(const Vec2Type& v) const
	{
		return checkCoords(static_cast<int32_t>(v.x), static_cast<int32_t>(v.y));
	}

	bool checkCoords(int32_t x, int32_t y) const
	{
		return static_cast<int32_t>(x) > 0 && static_cast<int32_t>(x) < (width - 1) &&
               static_cast<int32_t>(y) > 0 && static_cast<int32_t>(y) < (height - 1);
	}

	const T& get(int32_t x, int32_t y) const
	{
		return data[y * width + x];
	}

	template<typename Vec2Type>
	T& get(const Vec2Type& v)
	{
		return get(static_cast<int32_t>(v.x), static_cast<int32_t>(v.y));
	}

	template<typename Vec2Type>
	const T& get(const Vec2Type& v) const
	{
		return get(static_cast<int32_t>(v.x), static_cast<int32_t>(v.y));
	}

	template<typename Vec2Type>
	const T& getWrap(Vec2Type v) const
	{
		return getWrap(v.x, v.y);
	}

	const T& getWrap(int32_t x, int32_t y) const
	{
		return get(mod(x, width), mod(y, height));
	}

	T& get(int32_t x, int32_t y)
	{
		return data[y * width + x];
	}

	template<typename Vec2Type>
	void set(const Vec2Type& v, const T& obj)
	{
		set(v.x, v.y, obj);
	}

	void set(int32_t x, int32_t y, const T& obj)
	{
		data[y * width + x] = obj;
	}
};

struct CollisionCell
{
    static constexpr uint8_t cell_capacity = 4;
    static constexpr uint8_t max_cell_idx  = cell_capacity - 1;

    // Overlap workaround
	uint32_t objects_count              = 0;
    uint32_t objects[cell_capacity] = {};

	CollisionCell() = default;

	void addAtom(uint32_t id)
	{
        objects[objects_count] = id;
        objects_count += objects_count < max_cell_idx;
	}

	void clear()
	{
		objects_count = 0u;
	}

    void remove(uint32_t id)
    {
        for (uint32_t i{0}; i < objects_count; ++i) {
            if (objects[i] == id) {
                // Swap pop
                objects[i] = objects[objects_count - 1];
                --objects_count;
                return;
            }
        }
    }
};

struct CollisionGrid : public Grid<CollisionCell>
{
	CollisionGrid()
		: Grid<CollisionCell>()
	{}

	CollisionGrid(int32_t width, int32_t height)
		: Grid<CollisionCell>(width, height)
	{}

	bool addAtom(uint32_t x, uint32_t y, uint32_t atom)
	{
		const uint32_t id = x * height + y;
		// Add to grid
		data[id].addAtom(atom);
		return true;
	}

	void clear()
	{
		for (auto& c : data) {
            c.objects_count = 0;
        }
	}
};
