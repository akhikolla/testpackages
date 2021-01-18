/* Copyright 2012. Bloomberg Finance L.P.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:  The above
 * copyright notice and this permission notice shall be included in all copies
 * or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */
// blpapi_sessionoptions.h                                            -*-C++-*-
#ifndef INCLUDED_BLPAPI_SESSIONOPTIONS
#define INCLUDED_BLPAPI_SESSIONOPTIONS

//@PURPOSE: A common interface shared between publish and consumer sessions.
//
//@CLASSES:
//  blpapi::SessionOptions: user specified options when creating a session.
//
//@SEE_ALSO: blpapi_abstractsession, blpapi_session, blpapi_providersession
//
//@DESCRIPTION: This file defines a 'SessionOptions' class which is used
// to specify various options during session creation.
//
//
///Usage
///-----
// The following snippet shows to use the SessionOptions when creating a
// 'Session'.
//..
// #include <blpapi_session.h>
// SessionOptions sessionOptions;
// sessionOptions.setServerHost("127.0.0.1");
// Session session(sessionOptions);
// if (!session.start()) {
//      std::cout << "Failed to start session." << std::endl;
//      return;
// }
//..

#ifndef INCLUDED_BLPAPI_CALL
#include <blpapi_call.h>
#endif

#ifndef INCLUDED_BLPAPI_DEFS
#include <blpapi_defs.h>
#endif

#ifndef INCLUDED_BLPAPI_EXCEPTION
#include <blpapi_exception.h>
#endif

#ifndef INCLUDED_BLPAPI_TYPES
#include <blpapi_types.h>
#endif

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

BLPAPI_EXPORT
blpapi_SessionOptions_t *blpapi_SessionOptions_create(void);

BLPAPI_EXPORT
blpapi_SessionOptions_t *blpapi_SessionOptions_duplicate(
        const blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
void blpapi_SessionOptions_copy(blpapi_SessionOptions_t       *lhs,
                                const blpapi_SessionOptions_t *rhs);

BLPAPI_EXPORT
void blpapi_SessionOptions_destroy(blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_setServerHost(blpapi_SessionOptions_t *parameters,
                                        const char              *serverHost);

BLPAPI_EXPORT
int blpapi_SessionOptions_setServerPort(blpapi_SessionOptions_t *parameters,
                                        unsigned short           serverPort);

BLPAPI_EXPORT
int blpapi_SessionOptions_setServerAddress(blpapi_SessionOptions_t *parameters,
                                           const char              *serverHost,
                                           unsigned short           serverPort,
                                           size_t                   index);

BLPAPI_EXPORT
int blpapi_SessionOptions_removeServerAddress(
                                           blpapi_SessionOptions_t *parameters,
                                           size_t                   index);

BLPAPI_EXPORT
int blpapi_SessionOptions_setConnectTimeout(
                               blpapi_SessionOptions_t *parameters,
                               unsigned int             timeoutInMilliseconds);

BLPAPI_EXPORT
int blpapi_SessionOptions_setDefaultServices(
                                     blpapi_SessionOptions_t *parameters,
                                     const char              *defaultServices);

BLPAPI_EXPORT
int blpapi_SessionOptions_setDefaultSubscriptionService(
                                   blpapi_SessionOptions_t *parameters,
                                   const char              *serviceIdentifier);

BLPAPI_EXPORT
void blpapi_SessionOptions_setDefaultTopicPrefix(
                                           blpapi_SessionOptions_t *parameters,
                                           const char              *prefix);

BLPAPI_EXPORT
void blpapi_SessionOptions_setAllowMultipleCorrelatorsPerMsg(
                      blpapi_SessionOptions_t *parameters,
                      int                      allowMultipleCorrelatorsPerMsg);

BLPAPI_EXPORT
void blpapi_SessionOptions_setClientMode(blpapi_SessionOptions_t *parameters,
                                         int                      clientMode);

BLPAPI_EXPORT
void blpapi_SessionOptions_setMaxPendingRequests(
                                  blpapi_SessionOptions_t *parameters,
                                  int                      maxPendingRequests);

BLPAPI_EXPORT
void blpapi_SessionOptions_setAutoRestartOnDisconnection(
                                         blpapi_SessionOptions_t *parameters,
                                         int                      autoRestart);

BLPAPI_EXPORT
void blpapi_SessionOptions_setAutoRestart(
                                         blpapi_SessionOptions_t *parameters,
                                         int                      autoRestart);

BLPAPI_EXPORT
void blpapi_SessionOptions_setAuthenticationOptions(
                                         blpapi_SessionOptions_t *parameters,
                                         const char              *authOptions);

BLPAPI_EXPORT
void blpapi_SessionOptions_setNumStartAttempts(
                                    blpapi_SessionOptions_t *parameters,
                                    int                      numStartAttempts);

BLPAPI_EXPORT
void blpapi_SessionOptions_setMaxEventQueueSize(
                                   blpapi_SessionOptions_t *parameters,
                                   size_t                   maxEventQueueSize);

BLPAPI_EXPORT
int blpapi_SessionOptions_setSlowConsumerWarningHiWaterMark(
                                         blpapi_SessionOptions_t *parameters,
                                         float                    hiWaterMark);

BLPAPI_EXPORT
int blpapi_SessionOptions_setSlowConsumerWarningLoWaterMark(
                                         blpapi_SessionOptions_t *parameters,
                                         float                    loWaterMark);

BLPAPI_EXPORT
int blpapi_SessionOptions_setDefaultKeepAliveInactivityTime(
                                     blpapi_SessionOptions_t *parameters,
                                     int                      inactivityMsecs);

BLPAPI_EXPORT
int blpapi_SessionOptions_setDefaultKeepAliveResponseTimeout(
                                        blpapi_SessionOptions_t *parameters,
                                        int                      timeoutMsecs);

BLPAPI_EXPORT
int blpapi_SessionOptions_setKeepAliveEnabled(
                                           blpapi_SessionOptions_t *parameters,
                                           int                      isEnabled);

BLPAPI_EXPORT
void blpapi_SessionOptions_setRecordSubscriptionDataReceiveTimes(
                                        blpapi_SessionOptions_t *parameters,
                                        int                      shouldRecord);

BLPAPI_EXPORT
const char *blpapi_SessionOptions_serverHost(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
unsigned int blpapi_SessionOptions_serverPort(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_numServerAddresses(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_getServerAddress(
                                          blpapi_SessionOptions_t  *parameters,
                                          const char              **serverHost,
                                          unsigned short           *serverPort,
                                          size_t                    index);

BLPAPI_EXPORT
unsigned int blpapi_SessionOptions_connectTimeout(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
const char *blpapi_SessionOptions_defaultServices(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
const char *blpapi_SessionOptions_defaultSubscriptionService(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
const char *blpapi_SessionOptions_defaultTopicPrefix(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_allowMultipleCorrelatorsPerMsg(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_clientMode(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_maxPendingRequests(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_autoRestartOnDisconnection(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_autoRestart(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
const char *blpapi_SessionOptions_authenticationOptions(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_numStartAttempts(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
size_t blpapi_SessionOptions_maxEventQueueSize(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
float blpapi_SessionOptions_slowConsumerWarningHiWaterMark(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
float blpapi_SessionOptions_slowConsumerWarningLoWaterMark(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_defaultKeepAliveInactivityTime(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_defaultKeepAliveResponseTimeout(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_keepAliveEnabled(
                                          blpapi_SessionOptions_t *parameters);

BLPAPI_EXPORT
int blpapi_SessionOptions_recordSubscriptionDataReceiveTimes(
                                          blpapi_SessionOptions_t *parameters);

#ifdef __cplusplus
}

namespace BloombergLP {
namespace blpapi {
                         // ====================
                         // class SessionOptions
                         // ====================

class SessionOptions {
    // Contains the options which the user can specify when creating a
    // session.
    //
    // To use non-default options on a Session, create a
    // SessionOptions instance and set the required options and then
    // supply it when creating a Session.

    blpapi_SessionOptions_t *d_handle_p;

  public:

    // The possible options for how to connect to the API

    enum ClientMode {
        AUTO = BLPAPI_CLIENTMODE_AUTO,  // Automatic (desktop if available
                                        // otherwise server)
        DAPI = BLPAPI_CLIENTMODE_DAPI,  // Always connect to the desktop API
        SAPI = BLPAPI_CLIENTMODE_SAPI   // Always connect to the server API
    };

    SessionOptions();
        // Create a SessionOptions with all options set to the
        // standard defaults.

    SessionOptions(const SessionOptions& original);
        // Copy constructor

    ~SessionOptions();
        // Destroy this SessionOptions.

    // MANIPULATORS

    SessionOptions& operator=(const SessionOptions& rhs);
        // Assign to this object the value of the specified 'rhs' object.

    void setServerHost(const char *host);
        // Set the API server host to connect to when using the server
        // API to the specified 'host'.  A hostname or an IPv4 address
        // (that is, a.b.c.d).  The default is "127.0.0.1".

    void setServerPort(unsigned short port);
        // Set the port to connect to when using the server API to the
        // specified 'port'.  The default is 8194.

    int setServerAddress(
            const char     *serverHost,
            unsigned short  serverPort,
            size_t          index);
        // Set the server address at the specified 'index' using the specified
        // 'serverHost' and 'serverPort'.

    int removeServerAddress(size_t index);
        // Remove the server address at the specified 'index'.

    void setConnectTimeout(unsigned int timeoutMilliSeconds);
        // Set the connection timeout in milliseconds when connecting to the
        // API.  The default is 5000 milliseconds.  Behavior is not defined
        // unless the specified 'timeoutMilliSeconds' is in range of
        // [1 .. 120000] milliseconds

    void setDefaultServices(const char *defaultServices);
        // DEPRECATED
        // Set the default service for the session. This function is
        // deprecated; see 'setDefaultSubscriptionService'.

    void setDefaultSubscriptionService(const char *serviceIdentifier);
        // Set the default service for subscriptions which do not specify a
        // subscription server to the specified 'serviceIdentifier'. The
        // behavior is undefined unless 'serviceIdentifier' matches the regular
        // expression '^//[-_.a-zA-Z0-9]+/[-_.a-zA-Z0-9]+$'. The default is
        // "//blp/mktdata".  For more information on when this will be used see
        // 'QUALIFYING SUBSCRIPTION STRINGS' section in
        // 'blpapi_subscriptionlist'.

    void setDefaultTopicPrefix(const char *prefix);
        // Set the default topic prefix to be used when a subscription
        // does not specify a prefix to the specified 'prefix'. The default is
        // "/ticker/".  For more information on when this will be used see
        // 'QUALIFYING SUBSCRIPTION STRINGS' section in
        // 'blpapi_subscriptionlist'.

    void setAllowMultipleCorrelatorsPerMsg(
                                          bool allowMultipleCorrelatorsPerMsg);
        // Set whether the Session is allowed to associate more than
        // one CorrelationId with a Message to the specified
        // 'allowMultipleCorrelatorsPerMsg'.  The default is
        // 'false'.  This means that if you have multiple subscriptions
        // which overlap (that is a particular Message is relevant to
        // all of them) you will be presented with the same message
        // multiple times when you use the MessageIterator, each time
        // with a different CorrelationId.  If you specify 'true' for
        // this then a Message may be presented with multiple
        // CorrelationId's.

    void setClientMode(int clientMode);
        // Set how to connect to the API.  The default is AUTO which will try
        // to connect to the desktop API but fall back to the server API if the
        // desktop is not available.  DAPI always connects to the desktop API
        // and will fail if it is not available.  SAPI always connects to the
        // server API and will fail if it is not available.

    void setMaxPendingRequests(int maxPendingRequests);
        // Set the maximum number of requests which can be pending to
        // the specified 'maxPendingRequests'.  The default is 1024.

    void setAuthenticationOptions(const char *authOptions);
        // Set the specified 'authOptions' as authentication option.

    void setNumStartAttempts(int numStartAttempts);
        // Set the maximum number of attempts to start a session
        // by connecting a server.

    void setAutoRestartOnDisconnection(bool autoRestart);
        // Set whether automatically restarting connection if disconnected.

    void setMaxEventQueueSize(size_t eventQueueSize);
        // Set the maximum number of outstanding undelivered events per session
        // to the specified 'eventQueueSize'.  All subsequent events
        // delivered over the network will be dropped by the session if the
        // number of outstanding undelivered events is 'eventQueueSize',
        // the specified threshold.  The default value is 10000.

    void setSlowConsumerWarningHiWaterMark(float hiWaterMark);
        // Set the point at which "slow consumer" events will be generated,
        // using the specified 'highWaterMark' as a fraction of
        // 'maxEventQueueSize'; the default value is 0.75.  A warning event
        // will be generated when the number of outstanding undelivered events
        // passes above 'hiWaterMark * maxEventQueueSize'.  The behavior of the
        // function is undefined unless '0.0 < hiWaterMark <= 1.0'.  Further,
        // at the time that 'Session.start()' is called, it must be the case
        // that 'slowConsumerWarningLoWaterMark() * maxEventQueueSize()' <
        // 'slowConsumerWarningHiWaterMark() * maxEventQueueSize()'.

    void setSlowConsumerWarningLoWaterMark(float loWaterMark);
        // Set the point at which "slow consumer cleared" events will be
        // generated, using the specified 'loWaterMark' as a fraction of
        // 'maxEventQueueSize'; the default value is 0.5.  A warning cleared
        // event will be generated when the number of outstanding undelivered
        // events drops below 'loWaterMark * maxEventQueueSize'.
        // The behavior of the function is undefined unless
        // '0.0 <= loWaterMark < 1.0'.  Further, at the time that
        // 'Session.start()' is called, it must be the case that
        // 'slowConsumerWarningLoWaterMark() * maxEventQueueSize()' <
        // 'slowConsumerWarningHiWaterMark() * maxEventQueueSize()'.

    void setDefaultKeepAliveInactivityTime(int inactivityMsecs);
        // Set to the specified 'inactivityMsecs' the amount of time that no
        // traffic can be received on a connection before the ping-based
        // keep-alive mechanism is triggered; if no traffic is received for
        // this duration then a keep-alive ping is sent to the remote end to
        // solicit a response.  If 'inactivityMsecs == 0', then no keep-alive
        // pings will be sent.  The behavior of this function is undefined
        // unless 'inactivityMsecs' is a non-negative value.  The default value
        // is 20,000 milliseconds.  Note that not all back-end connections
        // provide ping-based keep-alives; this option is ignored by such
        // connections.

    void setDefaultKeepAliveResponseTimeout(int timeoutMsecs);
        // When a keep-alive ping is sent, wait for the specified
        // 'timeoutMsecs' to receive traffic (of any kind) before terminating
        // the connection due to inactivity.  If 'timeoutMsecs == 0', then
        // connections are never terminated due to the absence of traffic after
        // a keep-alive ping.  The behavior of this function is undefined
        // unless 'timeoutMsecs' is a non-negative value.  The default value is
        // 5,000 milliseconds.  Note that not all back-end connections provide
        // support for ping-based keep-alives; this option is ignored by such
        // connections.

    void setKeepAliveEnabled(bool isEnabled);
        // If the specified 'isEnabled' is 'false', then disable all keep-alive
        // mechanisms, both from the client to the server and from the server
        // to the client; otherwise enable keep-alive pings both from the
        // client to the server (as configured by
        // 'setDefaultKeepAliveInactivityTime' and
        // 'setDefaultKeepAliveResponseTimeout' if the connection supports
        // ping-based keep-alives), and from the server to the client as
        // specified by the server configuration.

    void setRecordSubscriptionDataReceiveTimes(bool shouldRecord);
        // Set whether the receipt time (accessed via
        // blpapi::Message::timeReceived) should be recorded for subscription
        // data messages. By default, the receipt time for these messages is
        // not recorded.

    // ACCESSORS
    const char *serverHost() const;
        // Return a pointer to the value of the server host option in this
        // SessionOptions instance.  The pointer is valid until this
        // SessionOptions is destroyed setServerHost() is called.

    unsigned short serverPort() const;
        // Return the server port that this session connects to.  If more than
        // one server addresses are specified, return the port of the first
        // server address.

    size_t numServerAddresses() const;
        // Return the number of server addresses.

    int getServerAddress(const char     **serverHost,
                         unsigned short  *serverPort,
                         size_t           index) const;
        // Put the server name and port into 'serverHost' and 'serverPort'
        // indexed by 'index'.  Return 0 if succeeded; otherwise, return
        // non-zero.

    unsigned int connectTimeout() const;
        // Return the value of the connection timeout option in this
        // SessionOptions instance in milliseconds.

    const char *defaultServices() const;
        // Return all default services in one string.

    const char *defaultSubscriptionService() const;
        // Return a pointer to the value of the default subscription
        // service option in this SessionOptions instance.  The pointer is
        // valid until this SessionOptions is destroyed or
        // setDefaultSubscriptionService() is called.

    const char *defaultTopicPrefix() const;
        // Return a pointer to the value of the default topic prefix
        // option in this SessionOptions instance.  The pointer is
        // valid until this SessionOptions is destroyed or
        // setDefaultTopicPrefix() is called.

    bool allowMultipleCorrelatorsPerMsg() const;
        // Return the value of the allow multiple correlators per
        // message option in this SessionOptions instance.

    int clientMode() const;
        // Return the value of the client mode option in this
        // SessionOptions instance.

    int maxPendingRequests() const;
        // Return the value of the maximum pending request option in
        // this SessionOptions instance.

    blpapi_SessionOptions_t *handle() const;
        // Return the handle of the current session

    bool autoRestartOnDisconnection() const;
        // Return whether automatically restarting connection if disconnected.

    const char *authenticationOptions() const;
        // Return authentication options in a string.

    int numStartAttempts() const;
        // Return the maximum number of attempts to start a session
        // by connecting a server.

    size_t maxEventQueueSize() const;
        // Return the value of maximum outstanding undelivered events
        // that the session is configured with.

    float slowConsumerWarningHiWaterMark() const;
        // Return the fraction of maxEventQueueSize at which "slow consumer"
        // event will be generated.

    float slowConsumerWarningLoWaterMark() const;
        // Return the fraction of maxEventQueueSize at which
        // "slow consumer cleared" event will be generated.

    int defaultKeepAliveInactivityTime() const;
        // Return the interval (in milliseconds) a connection has to remain
        // inactive (receive no data) before a keep alive probe will be sent.

    int defaultKeepAliveResponseTimeout() const;
        // Return the time (in milliseconds) the library will wait for response
        // to a keep alive probe before declaring it lost.

    bool keepAliveEnabled() const;
        // Return 'true' if the keep-alive mechanism is enabled; otherwise
        // return 'false'.

    bool recordSubscriptionDataReceiveTimes() const;
        // Return whether the receipt time (accessed via
        // blpapi::Message::timeReceived) should be recorded for subscription
        // data messages.
};

// ============================================================================
//                      INLINE FUNCTION DEFINITIONS
// ============================================================================

                            // --------------------
                            // class SessionOptions
                            // --------------------
inline
SessionOptions::SessionOptions()
{
    d_handle_p = blpapi_SessionOptions_create();
#if BLPAPI_COMPAT_33X
    blpapi_SessionOptions_setClientMode(
            d_handle_p,
            BLPAPI_CLIENTMODE_AUTO | BLPAPI_CLIENTMODE_COMPAT_33X);
#endif
}

inline
SessionOptions::SessionOptions(const SessionOptions& options)
{
    d_handle_p = blpapi_SessionOptions_duplicate(options.handle());
}

inline
SessionOptions::~SessionOptions()
{
    blpapi_SessionOptions_destroy(d_handle_p);
}

inline
SessionOptions& SessionOptions::operator=(const SessionOptions& rhs)
{
    blpapi_SessionOptions_copy(this->handle(), rhs.handle());
    return *this;
}

inline
void SessionOptions::setServerHost(const char *newServerHost)
{
    blpapi_SessionOptions_setServerHost(d_handle_p, newServerHost);
}

inline
void SessionOptions::setServerPort(unsigned short newServerPort)
{
    blpapi_SessionOptions_setServerPort(d_handle_p, newServerPort);
}

inline
int SessionOptions::setServerAddress(const char     *newServerHost,
                                     unsigned short  newServerPort,
                                     size_t          index)
{
    return blpapi_SessionOptions_setServerAddress(d_handle_p,
                                                  newServerHost,
                                                  newServerPort,
                                                  index);
}

inline
int SessionOptions::removeServerAddress(size_t index)
{
    return blpapi_SessionOptions_removeServerAddress(d_handle_p, index);
}

inline
void SessionOptions::setConnectTimeout(unsigned int timeoutMilliSeconds)
{
    ExceptionUtil::throwOnError(
        blpapi_SessionOptions_setConnectTimeout(
            d_handle_p, timeoutMilliSeconds));
}

inline
void SessionOptions::setDefaultServices(const char *newDefaultServices)
{
    blpapi_SessionOptions_setDefaultServices(d_handle_p, newDefaultServices);
}

inline
void SessionOptions::setDefaultSubscriptionService(
                                                 const char *serviceIdentifier)
{
    blpapi_SessionOptions_setDefaultSubscriptionService(
            d_handle_p,
            serviceIdentifier);
}

inline
void SessionOptions::setDefaultTopicPrefix(const char *prefix)
{
    blpapi_SessionOptions_setDefaultTopicPrefix(
            d_handle_p,
            prefix);
}

inline
void SessionOptions::setAllowMultipleCorrelatorsPerMsg(
                                        bool newAllowMultipleCorrelatorsPerMsg)
{
    blpapi_SessionOptions_setAllowMultipleCorrelatorsPerMsg(
            d_handle_p,
            newAllowMultipleCorrelatorsPerMsg);
}

inline
void SessionOptions::setClientMode(int newClientMode)
{
#if BLPAPI_COMPAT_33X
    newClientMode |= BLPAPI_CLIENTMODE_COMPAT_33X;
#endif

    blpapi_SessionOptions_setClientMode(
            d_handle_p,
            newClientMode);
}

inline
void SessionOptions::setMaxPendingRequests(int newMaxPendingRequests)
{
    blpapi_SessionOptions_setMaxPendingRequests(
            d_handle_p,
            newMaxPendingRequests);
}

inline
void SessionOptions::setAutoRestartOnDisconnection(bool autoRestart)
{
    blpapi_SessionOptions_setAutoRestartOnDisconnection(
            d_handle_p,
            autoRestart? 1: 0);
}

inline
void SessionOptions::setAuthenticationOptions(const char *authOptions)
{
    blpapi_SessionOptions_setAuthenticationOptions(
            d_handle_p, authOptions);
}

inline
void SessionOptions::setNumStartAttempts(int newNumStartAttempts)
{
    blpapi_SessionOptions_setNumStartAttempts(
            d_handle_p, newNumStartAttempts);
}

inline
void SessionOptions::setMaxEventQueueSize(size_t eventQueueSize)
{
    BLPAPI_CALL_SESSIONOPTIONS_SETMAXEVENTQUEUESIZE(
            d_handle_p,
            eventQueueSize);
}

inline
void SessionOptions::setSlowConsumerWarningHiWaterMark(float hiWaterMark)
{
    ExceptionUtil::throwOnError(
        BLPAPI_CALL_SESSIONOPTIONS_SETSLOWCONSUMERHIWATERMARK(
            d_handle_p,
            hiWaterMark));
}

inline
void SessionOptions::setSlowConsumerWarningLoWaterMark(float loWaterMark)
{
    ExceptionUtil::throwOnError(
        BLPAPI_CALL_SESSIONOPTIONS_SETSLOWCONSUMERLOWATERMARK(
            d_handle_p,
            loWaterMark));
}

inline
void SessionOptions::setDefaultKeepAliveInactivityTime(int inactivityTime)
{
    ExceptionUtil::throwOnError(
        BLPAPI_CALL_SESSIONOPTIONS_SETDEFAULTKEEPALIVEINACTIVITYTIME(
            d_handle_p,
            inactivityTime));
}

inline
void SessionOptions::setDefaultKeepAliveResponseTimeout(int responseTimeout)
{
    ExceptionUtil::throwOnError(
        BLPAPI_CALL_SESSIONOPTIONS_SETDEFAULTKEEPALIVERESPONSETIMEOUT(
            d_handle_p,
            responseTimeout));
}

inline
void SessionOptions::setKeepAliveEnabled(bool isEnabled)
{
    ExceptionUtil::throwOnError(
        BLPAPI_CALL_SESSIONOPTIONS_SETKEEPALIVEENABLED(d_handle_p, isEnabled));
}

inline
void SessionOptions::setRecordSubscriptionDataReceiveTimes(bool shouldRecrod)
{
    BLPAPI_CALL_SESSIONOPTION_SETRECORDSUBSCRIPTIONDATARECEIVETIMES(
            d_handle_p,
            shouldRecrod);
}

inline
const char *SessionOptions::serverHost() const
{
    return blpapi_SessionOptions_serverHost(d_handle_p);
}

inline
unsigned short SessionOptions::serverPort() const
{
    return static_cast<unsigned short>(
        blpapi_SessionOptions_serverPort(d_handle_p));
}

inline
size_t SessionOptions::numServerAddresses() const
{
    return blpapi_SessionOptions_numServerAddresses(d_handle_p);
}

inline
int SessionOptions::getServerAddress(
                            const char     **serverHostOut,
                            unsigned short  *serverPortOut,
                            size_t           index) const
{
    return blpapi_SessionOptions_getServerAddress(d_handle_p,
                                                  serverHostOut,
                                                  serverPortOut,
                                                  index);
}

inline
unsigned int SessionOptions::connectTimeout() const
{
    return blpapi_SessionOptions_connectTimeout(d_handle_p);
}

inline
const char *SessionOptions::defaultServices() const
{
    return blpapi_SessionOptions_defaultServices(d_handle_p);
}

inline
const char *SessionOptions::defaultSubscriptionService() const
{
    return blpapi_SessionOptions_defaultSubscriptionService(d_handle_p);
}

inline
const char *SessionOptions::defaultTopicPrefix() const
{
    return blpapi_SessionOptions_defaultTopicPrefix(d_handle_p);
}

inline
bool SessionOptions::allowMultipleCorrelatorsPerMsg() const
{
    return blpapi_SessionOptions_allowMultipleCorrelatorsPerMsg(d_handle_p) ?
        true : false;
}

inline
int SessionOptions::clientMode() const
{
    return blpapi_SessionOptions_clientMode(d_handle_p);
}

inline
int SessionOptions::maxPendingRequests() const
{
    return blpapi_SessionOptions_maxPendingRequests(d_handle_p);
}

inline
bool SessionOptions::autoRestartOnDisconnection() const
{
    return blpapi_SessionOptions_autoRestartOnDisconnection(d_handle_p) != 0;
}

inline
const char *SessionOptions::authenticationOptions() const
{
    return blpapi_SessionOptions_authenticationOptions(d_handle_p);
}

inline
int SessionOptions::numStartAttempts() const
{
    return blpapi_SessionOptions_numStartAttempts(d_handle_p);
}

inline
size_t SessionOptions::maxEventQueueSize() const
{
    return BLPAPI_CALL_SESSIONOPTIONS_MAXEVENTQUEUESIZE(d_handle_p);
}

inline
float SessionOptions::slowConsumerWarningHiWaterMark() const
{
    return BLPAPI_CALL_SESSIONOPTIONS_SLOWCONSUMERHIWATERMARK(d_handle_p);
}

inline
float SessionOptions::slowConsumerWarningLoWaterMark() const
{
    return BLPAPI_CALL_SESSIONOPTIONS_SLOWCONSUMERLOWATERMARK(d_handle_p);
}

inline
int SessionOptions::defaultKeepAliveInactivityTime() const
{
    return
        BLPAPI_CALL_SESSIONOPTIONS_DEFAULTKEEPALIVEINACTIVITYTIME(d_handle_p);
}

inline
int SessionOptions::defaultKeepAliveResponseTimeout() const
{
    return
        BLPAPI_CALL_SESSIONOPTIONS_DEFAULTKEEPALIVERESPONSETIMEOUT(d_handle_p);
}

inline
bool SessionOptions::keepAliveEnabled() const
{
    return BLPAPI_CALL_SESSIONOPTIONS_KEEPALIVEENABLED(d_handle_p) != 0
        ? true
        : false;
}

inline
bool SessionOptions::recordSubscriptionDataReceiveTimes() const
{
    return
        BLPAPI_CALL_SESSIONOPTION_RECORDSUBSCRIPTIONDATARECEIVETIMES(
                                                                    d_handle_p)
        ? true
        : false;
}

inline
blpapi_SessionOptions_t *SessionOptions::handle() const
{
    return d_handle_p;
}

}  // close namespace blpapi
}  // close namespace BloombergLP

#endif // #ifdef __cplusplus
#endif // #ifndef INCLUDED_BLPAPI_SESSIONOPTIONS
