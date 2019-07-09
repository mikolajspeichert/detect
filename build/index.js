if (process.env.NODE_ENV === 'production') {
  module.exports = require('./detect.production.js')
} else {
  module.exports = require('./detect.development.js')
}
