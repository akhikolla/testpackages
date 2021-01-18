#ifndef ROBOTS_CPP_H
#define ROBOTS_CPP_H

#include <sstream>
#include <unordered_map>
#include <vector>

#include "agent.h"

namespace Rep
{

    class Robots
    {
    public:
        typedef std::unordered_map<std::string, Agent> agent_map_t;
        typedef std::vector<std::string> sitemaps_t;
        agent_map_t agents_;

        /**
         * Create a robots.txt from a utf-8-encoded string.
         */
        Robots(const std::string& content);

        /**
         * Instantiate a Robots object.
         */
        Robots(
            const agent_map_t& agents,
            const sitemaps_t& sitemaps)
            : agents_(agents)
            , sitemaps_(sitemaps)
            , default_(agents_["*"]) {}

        /**
         * Get the sitemaps in this robots.txt
         */
        const sitemaps_t& sitemaps() const { return sitemaps_; }

        /**
         * Get the agent with the corresponding name.
         */
        const Agent& agent(const std::string& name) const;

        /**
         * Return true if agent is allowed to fetch the URL (either a
         * full URL or a path).
         */
        bool allowed(const std::string& path, const std::string& name) const;

        std::string str() const;

        /**
         * Return the robots.txt URL corresponding to the provided URL.
         */
        static std::string robotsUrl(const std::string& url);

    private:
        static void strip(std::string& string);

        static bool getpair(
            std::istringstream& stream, std::string& key, std::string& value);

        sitemaps_t sitemaps_;
        Agent& default_;
    };
}

#endif
